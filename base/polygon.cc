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

/*!\file polygons.cc
   \brief Generate polyhedra based on polygons.
*/

#include <algorithm>
#include <float.h>
#include <stdlib.h>
#include <vector>

#include "polygon.h"
#include "private_std_polys.h"
#include "symmetry.h"
#include "utils.h"

using std::vector;

namespace anti {

Polygon::Polygon(int N, int D, int type, int subtype)
    : radius(vector<double>(2, NAN)), height(vector<double>(2, NAN)),
      edge(vector<double>(2, NAN)), twist_angle(vector<double>(2, NAN)),
      type(type), subtype(subtype)
{
  set_fraction(N, D);
}

Status Polygon::set_type(int typ)
{
  if (type >= unknown && type < types_end) {
    type = typ;
    subtype = sub_default;
    return Status::ok();
  }
  else
    return Status::error("unknown model type");
}

Status Polygon::set_subtype(int subtyp)
{
  unsigned int params;
  vector<unsigned int> conflicts;
  int max_subtype = get_acceptable_params(params, conflicts);

  if (subtyp >= 0 && subtyp <= max_subtype) {
    subtype = subtyp;
    return Status::ok();
  }
  else
    return Status::error("subtype not valid for this model type");
}

Status Polygon::set_fraction(int N, int D)
{
  if (abs(N) < 2)
    return Status::error("number of sides must be an integer 2 or greater");

  if (D < 1)
    return Status::error(
        "denominator of polygon fraction must be 1 or greater");

  if (D % N == 0)
    return Status::error("denominator of polygon fraction cannot be a "
                         "multiple of the number of sides");

  parts = gcd(N, D);
  num_sides = N / parts;
  step = D / parts;

  return Status::ok();
}

int Polygon::get_acceptable_params(unsigned int &params,
                                   vector<unsigned int> &conflicts)
{
  int max_subtype = 0;
  params = R0 | H0 | E0;
  conflicts.clear();
  conflicts.push_back(R0 | E0);
  conflicts.push_back(H0 | E1);
  enum {
    unknown = 0,
    prism,
    antiprism,
    pyramid,
    dipyramid,
    cupola,
    orthobicupola,
    gyrobicupola,
    snub_antiprism,
    dihedron
  };
  params = R0 | H0 | E0;
  switch (type) {
  case prism:
    max_subtype = 3;
    params |= A0;
    if (subtype == sub_prism_crown)
      params |= A1;
    break;
  case antiprism:
    max_subtype = 5;
    params |= A0;
    if (subtype == sub_antiprism_scalenohedron ||
        subtype == sub_antiprism_subdivided_scalenohedron)
      params |= E1;
    if (subtype == sub_antiprism_crown)
      params |= A1;
    break;
  case pyramid:
    max_subtype = 3;
    params |= A0;
    if (subtype == sub_default || subtype == sub_pyramid_elongated ||
        subtype == sub_pyramid_gyroelongated)
      params |= E1;
    break;
  case dipyramid:
    max_subtype = 4;
    params |= A0;
    if (subtype == sub_default || subtype == sub_dipyramid_elongated ||
        subtype == sub_dipyramid_gyroelongated)
      params |= E1;
    else if (subtype == sub_dipyramid_scalenohedron)
      params |= R1;
    break;
  case cupola:
    max_subtype = 3;
    if (subtype == sub_cupola_elongated || subtype == sub_cupola_gyroelongated)
      params |= E1;
    break;
  case orthobicupola:
    max_subtype = 2;
    if (subtype == sub_cupola_elongated || subtype == sub_cupola_gyroelongated)
      params |= E1;
    break;
  case gyrobicupola:
    max_subtype = 2;
    if (subtype == sub_cupola_elongated || subtype == sub_cupola_gyroelongated)
      params |= E1;
    break;
  case dihedron:
    max_subtype = 1;
    conflicts.clear();
    break;
  case snub_antiprism:
    max_subtype = 1;
    break;
  case crown:
    params |= A0 | A1 | E1;
    conflicts.push_back(A0 | A1);
    max_subtype = 0;
    break;
  default:
    params = 0;
    conflicts.clear();
    max_subtype = -1;
    break;
  }

  return max_subtype; // >=0 if valid type
}

unsigned int Polygon::get_params_set()
{
  return R0 * value_is_set(radius[0]) | R1 * value_is_set(radius[1]) |
         H0 * value_is_set(height[0]) | H1 * value_is_set(height[1]) |
         E0 * value_is_set(edge[0]) | E1 * value_is_set(edge[1]) |
         A0 * value_is_set(twist_angle[0]) | A1 * value_is_set(twist_angle[1]);
}

bool Polygon::has_params(unsigned int params)
{
  if (!type)
    return false;

  unsigned int acceptable_params;
  vector<unsigned int> conflicts;
  get_acceptable_params(acceptable_params, conflicts);
  return acceptable_params & params;
}

double Polygon::get_polygon_radius()
{
  // prefer radius over edge
  if (value_is_set(radius[0]))
    return radius[0];
  else if (value_is_set(edge[0]))
    return edge[0] / (2 * sin(angle() / 2));
  else
    return 1.0;
}

double Polygon::get_polygon_edge()
{
  return get_polygon_radius() * (2 * sin(angle() / 2));
}

void Polygon::add_polygon(Geometry &geom, double ht)
{
  double rad = get_polygon_radius();
  vector<int> face(num_sides);
  int offset = geom.verts().size();
  for (int i = 0; i < num_sides; i++) {
    geom.add_vert(rad * cos(i * angle()) * Vec3d::X + ht * Vec3d::Z +
                  -rad * sin(i * angle()) * Vec3d::Y);
    face[i] = offset + i;
  }
  geom.add_face(face);
}

Status Polygon::make_poly(Geometry &geom)
{
  Geometry part;
  Status stat;
  if (type == prism)
    stat = make_prism_part(part);
  else if (type == antiprism)
    stat = make_antiprism_part(part);
  else if (type == pyramid)
    stat = make_pyramid_part(part);
  else if (type == dipyramid)
    stat = make_dipyramid_part(part);
  else if (type == cupola)
    stat = make_cupola_part(part);
  else if (type == orthobicupola)
    stat = make_orthobicupola_part(part);
  else if (type == gyrobicupola)
    stat = make_gyrobicupola_part(part);
  else if (type == dihedron)
    stat = make_dihedron_part(part);
  else if (type == snub_antiprism)
    stat = make_snub_antiprism_part(part);
  else if (type == crown) {
    if (value_is_set(get_twist_angle(1)))
      return make_crown_full(geom); // no repetitions needed
    stat = make_crown_part(part);
  }
  else
    stat.set_error(msg_str("unknown type (%d) for polygon-based model", type));

  if (stat) {
    part.orient();
    repeat_part(geom, part);
  }

  return stat;
}

void Polygon::repeat_part(Geometry &geom, const Geometry &part)
{
  geom.append(part);
  for (int i = 1; i < parts; i++) {
    Geometry rep = part;
    rep.transform(Trans3d::rotate(Vec3d::Z, 2 * M_PI * i / parts / num_sides));
    geom.append(rep);
  }
}

void Polygon::dump()
{
  unsigned int params;
  vector<unsigned int> conflicts;
  int max_subtype = get_acceptable_params(params, conflicts);

  fprintf(stderr, "\npolygon %d x {%d/%d}\n", parts, num_sides, step);
  fprintf(stderr, "type, %d, subtype %d (out of %d)\n", type, subtype,
          max_subtype);
  fprintf(stderr, "R=[%g, %g], H=[%g, %g], E=[%g, %g], A=[%g, %g]\n\n",
          radius[0], radius[1], height[0], height[1], edge[0], edge[1],
          twist_angle[0], twist_angle[1]);
}

//--------------------------------------------------------------------

Status Polygon::make_dihedron_part(Geometry &geom)
{
  enum { subtype_polygon = 1 };
  add_polygon(geom);
  if (subtype == sub_default) {
    geom.add_face(geom.faces(0));
    reverse(geom.raw_faces()[1].begin(), geom.raw_faces()[1].end());
  }

  return Status::ok();
}

Status Polygon::make_prism_crown_part(Geometry &geom)
{
  Polygon pri(*this);
  pri.set_subtype(0);
  pri.set_twist_angle(0, NAN);
  Status stat = pri.make_prism_part(geom);
  if (stat.is_error())
    return stat;
  geom.clear(FACES);
  int step2 = value_is_set(get_twist_angle(1))
                  ? 1
                  : int(floor(rad2deg(get_twist_angle(1)) + 0.5));
  step2 = ((step2 % num_sides) + num_sides) % num_sides;

  int N = num_sides;
  for (int i = 0; i < N; i++) {
    geom.add_face((i - step2 + num_sides) % num_sides,
                  (i + 1) % num_sides + num_sides, (i + 1 + step2) % num_sides,
                  (i) % num_sides + num_sides, -1);

    geom.add_face((i - step2 + num_sides) % num_sides + num_sides,
                  (i + 1) % num_sides, (i + 1 + step2) % num_sides + num_sides,
                  (i) % num_sides, -1);
  }

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_prism_part(Geometry &geom)
{
  double twist_ang = twist_angle[0];
  if (value_is_set(twist_ang) || subtype == sub_prism_trapezohedron ||
      subtype == sub_prism_antiprism) {
    Polygon ant(*this);
    ant.set_type(antiprism);
    if (value_is_set(twist_angle[0]) || subtype == sub_prism_antiprism)
      ant.set_subtype(sub_default);
    if (subtype == sub_prism_trapezohedron) // CHECK if ang but not trap!
      ant.set_subtype(sub_antiprism_trapezohedron);

    double ang = value_is_set(twist_ang) ? twist_ang : 0.0;
    ant.set_twist_angle(0, ang - angle() / 2);
    return ant.make_poly(geom);
  }
  else if (subtype == sub_prism_crown)
    return make_prism_crown_part(geom);

  double ht;
  if (value_is_set(height[0]))
    ht = height[0];
  else if (value_is_set(edge[1]))
    ht = edge[1];
  else
    ht = get_polygon_radius();

  vector<Vec3d> verts;
  vector<vector<int>> faces;
  verts.resize(2 * num_sides);
  faces.resize(2 + num_sides);
  Geometry pgon;
  add_polygon(pgon);
  for (int i = 0; i < num_sides; i++) {
    verts[i] = pgon.verts(i) + (ht / 2) * Vec3d::Z;
    verts[i + num_sides] = pgon.verts(i) - (ht / 2) * Vec3d::Z;
    faces[0].push_back(i);
    faces[1].push_back(i + num_sides);
    faces[2 + i].push_back(i);
    faces[2 + i].push_back((i + 1) % num_sides);
    faces[2 + i].push_back((i + 1) % num_sides + num_sides);
    faces[2 + i].push_back(i + num_sides);
  }
  reverse(faces[0].begin(), faces[0].end());
  geom.add_verts(verts);
  geom.add_faces(faces);

  return Status::ok();
}

//--------------------------------------------------------------------

double Polygon::get_antiprism_height()
{
  if (value_is_set(height[0]))
    return height[0];
  else if (value_is_set(edge[1])) {
    double dist = 2 * get_polygon_radius() * sin(angle() / 4);
    return (fabs(edge[1]) > fabs(dist)) ? sqrt(edge[1] * edge[1] - dist * dist)
                                        : 0.0;
  }
  else
    return get_polygon_radius();
}

Status Polygon::make_antiprism_scal_part(Geometry &geom)
{
  Polygon ant(*this);
  ant.set_type(antiprism);
  ant.set_subtype(sub_default);
  ant.set_twist_angle(0, 0.0);
  Status stat = ant.make_antiprism_part(geom);
  if (stat.is_error())
    return stat;
  geom.clear(FACES);

  double ht = get_antiprism_height();
  double apex_ht = value_is_set(height[1]) ? height[1] : 2 * ht;
  int apex_idx1 = geom.add_vert(Vec3d(0, 0, apex_ht));
  int apex_idx2 = geom.add_vert(Vec3d(0, 0, -apex_ht));

  for (int i = 0; i < num_sides; i++) {
    geom.add_face(i, i + num_sides, apex_idx1, -1);
    geom.add_face(i + num_sides, i, apex_idx2, -1);
    geom.add_face(i + num_sides, (i + 1) % num_sides, apex_idx1, -1);
    geom.add_face((i + 1) % num_sides, i + num_sides, apex_idx2, -1);
  }

  return Status::ok();
}

static void add_scal_faces(Geometry &geom, int v0, int v1, int v2, int v3,
                           double ht2)
{
  Vec3d cent = 0.25 * (geom.verts(v0) + geom.verts(v1) + geom.verts(v2) +
                       geom.verts(v3));
  Vec3d norm =
      vcross(geom.verts(v0) - geom.verts(v2), geom.verts(v3) - geom.verts(v1))
          .unit();
  int ap = geom.add_vert(cent + norm * ht2);
  geom.add_face(v0, v1, ap, -1);
  geom.add_face(v1, v2, ap, -1);
  geom.add_face(v2, v3, ap, -1);
  geom.add_face(v3, v0, ap, -1);
}

Status Polygon::make_antiprism_subscal_part(Geometry &geom)
{
  double ht = get_antiprism_height();
  double ht2 = (value_is_set(height[1])) ? height[1] : 0.0;
  double rad = get_polygon_radius();

  double E = sqrt(ht * ht + pow(2 * rad * sin(angle() / 4), 2));
  if (fabs(rad) > E - epsilon)
    return Status::error("antprism slant height to short to close "
                         "polyhedron at apex");

  double apex_ht = ht / 2; // default to "flat" on failure
  if (rad < E)
    apex_ht += sqrt(E * E - rad * rad);

  Polygon ant(*this);
  ant.set_type(antiprism);
  ant.set_subtype(sub_default);
  ant.set_twist_angle(0, 0.0);
  Status stat = ant.make_antiprism_part(geom);
  if (stat.is_error())
    return stat;
  geom.clear(FACES);

  int apex_idx1 = geom.add_vert(Vec3d(0, 0, apex_ht));
  int apex_idx2 = geom.add_vert(Vec3d(0, 0, -apex_ht));

  for (int i = 0; i < num_sides; i++) {
    add_scal_faces(geom, i, i + num_sides, (i + 1) % num_sides, apex_idx1, ht2);
    add_scal_faces(geom, (i + 1) % num_sides + num_sides, (i + 1) % num_sides,
                   i + num_sides, apex_idx2, ht2);
  }
  return Status::ok();
}

Status Polygon::make_antiprism_crown_part(Geometry &geom)
{
  Polygon ant(*this);
  ant.set_subtype(0);
  ant.set_twist_angle(0, 0.0);
  Status stat = ant.make_antiprism_part(geom);
  if (stat.is_error())
    return stat;
  geom.clear(FACES);
  int step2 = value_is_set(get_twist_angle(1))
                  ? 1
                  : int(floor(rad2deg(get_twist_angle(1)) + 0.5));
  step2 = ((step2 % num_sides) + num_sides) % num_sides;

  for (int i = 0; i < num_sides; i++) {
    geom.add_face((i + 1 - step2 + num_sides) % num_sides,
                  (i + 1) % num_sides + num_sides, (i + 1 + step2) % num_sides,
                  (i) % num_sides + num_sides, -1);
    geom.add_face((i + 1 - step2 + num_sides) % num_sides + num_sides,
                  (i + 1) % num_sides, (i + 1 + step2) % num_sides + num_sides,
                  (i + 2) % num_sides, -1);
  }

  return Status::ok();
}

Status Polygon::make_antiprism_part(Geometry &geom)
{
  if (subtype == sub_antiprism_scalenohedron)
    return make_antiprism_scal_part(geom);
  else if (subtype == sub_antiprism_subdivided_scalenohedron)
    return make_antiprism_subscal_part(geom);

  if (value_is_set(edge[1])) {
    double dist = 2 * get_polygon_radius() * sin(angle() / 4);
    if (fabs(edge[1]) - fabs(dist) < -epsilon)
      return Status::error("slant edge too short to reach between vertices");
  }
  else if (subtype == sub_antiprism_crown)
    return make_antiprism_crown_part(geom);

  double ht = get_antiprism_height();
  int extra_verts = 2 * (subtype == sub_antiprism_trapezohedron) +
                    1 * (subtype == sub_antiprism_antihermaphrodite);
  int extra_faces = 2 - extra_verts;

  vector<Vec3d> verts;
  vector<vector<int>> faces;
  vector<vector<int>> caps(2);
  verts.resize(2 * num_sides + extra_verts);
  faces.resize(2 * num_sides);

  Geometry pgon, pgon2;
  add_polygon(pgon);
  double twist_ang = value_is_set(twist_angle[0]) ? twist_angle[0] : 0.0;
  pgon.transform(Trans3d::rotate(Vec3d::Z, angle() / 4 - twist_ang / 2));
  add_polygon(pgon2);
  pgon2.transform(Trans3d::rotate(Vec3d::Z, -angle() / 4 + twist_ang / 2));
  if (extra_verts) {
    double apex_ht = 0.5 * ht;
    if (num_sides != 2)
      apex_ht *= (cos(angle() / 2) + cos(twist_ang)) /
                 (cos(twist_ang) - cos(angle() / 2));
    if (apex_ht > FLT_MAX)
      apex_ht = FLT_MAX;
    else if (apex_ht < -FLT_MAX)
      apex_ht = -FLT_MAX;

    verts[2 * num_sides] = Vec3d(0, 0, -apex_ht);
    if (extra_verts == 2)
      verts[2 * num_sides + 1] = Vec3d(0, 0, apex_ht);
  }

  for (int i = 0; i < num_sides; i++) {
    verts[i] = pgon.verts(i) + (ht / 2) * Vec3d::Z;
    verts[i + num_sides] = pgon2.verts(i) - (ht / 2) * Vec3d::Z;
    if (extra_faces) {
      caps[0].push_back(i);
      if (extra_faces > 1)
        caps[1].push_back(i + num_sides);
    }
    faces[i].push_back(i);
    if (extra_verts > 1)
      faces[i].push_back(2 * num_sides + 1);
    faces[i].push_back((i + 1) % num_sides);
    faces[i].push_back(i + num_sides);
    faces[i + num_sides].push_back((i + 1) % num_sides + num_sides);
    if (extra_verts)
      faces[i + num_sides].push_back(2 * num_sides);
    faces[i + num_sides].push_back(i % num_sides + num_sides);
    faces[i + num_sides].push_back((i + 1) % num_sides);
  }

  if (extra_faces) {
    reverse(caps[0].begin(), caps[0].end());
    geom.add_face(caps[0]);
    if (extra_faces > 1) {
      geom.add_face(caps[1]);
    }
  }
  geom.add_verts(verts);
  geom.add_faces(faces);

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_pyramid_part(Geometry &geom)
{
  double rad = get_polygon_radius();
  double ht;
  if (value_is_set(height[0]))
    ht = height[0];
  else if (value_is_set(edge[1])) {
    if (edge[1] - rad < -epsilon)
      return Status::error("slant edge too short to reach apex");

    ht = (edge[1] > rad) ? sqrt(edge[1] * edge[1] - rad * rad) : 0;
  }
  else
    ht = rad;

  if ((subtype == sub_default && value_is_set(twist_angle[0])) ||
      subtype == sub_pyramid_antihermaphrodite) {
    if (num_sides == 2)
      return Status::error(
          "cannot make an antihermaphrodite from a digonal pyramid");

    Polygon ant(*this);
    ant.set_type(antiprism);
    ant.set_subtype(sub_antiprism_antihermaphrodite);

    double twist_ang =
        value_is_set(twist_angle[0]) ? twist_angle[0] - angle() / 2 : 0.0;
    ant.set_twist_angle(0, twist_ang);
    double ant_ht = ht * (cos(twist_ang) / cos(angle() / 2) - 1);
    ant.set_height(0, ant_ht);
    ant.make_poly(geom);
    geom.verts(geom.verts().size() - 1) = Vec3d(0, 0, -(0.5 * ant_ht + ht));
    return Status::ok();
  }

  Geometry pyr;
  add_polygon(pyr);
  for (int i = 0; i < num_sides; i++) {
    vector<int> face(3);
    face[0] = i;
    face[1] = (i + 1) % num_sides;
    face[2] = num_sides;
    pyr.add_face(face);
  }
  pyr.add_vert(ht * Vec3d::Z);

  if (subtype == sub_default) {
    geom.append(pyr);
  }
  else if (subtype == sub_pyramid_elongated) {
    Polygon pri(*this);
    pri.set_type(prism);
    pri.set_subtype(sub_default);
    pri.set_height(0, value_is_set(height[1]) ? height[1] : get_polygon_edge());
    Status stat = pri.make_prism_part(geom);
    if (stat.is_error())
      return stat;
    if (num_sides > 2)
      face_bond(&geom, &pyr);
    else {
      geom.del(FACES, vector<int>(1, 0));
      int apex = geom.add_vert(Vec3d(0, 0, pri.get_height(0) / 2 + ht));
      geom.add_face(0, 1, apex, -1);
      geom.add_face(1, 0, apex, -1);
    }
  }
  else if (subtype == sub_pyramid_gyroelongated) {
    Polygon ant(*this);
    ant.set_type(antiprism);
    ant.set_subtype(sub_default);
    if (value_is_set(height[1]))
      ant.set_height(0, value_is_set(height[1])); // prefer second height
    else
      ant.set_edge(1, get_polygon_edge()); // else try polygon edge

    if (!ant.make_antiprism_part(geom)) { // may fail for pgon edge
      ant.set_height(0, ht);              // so use pyr height
      ant.make_antiprism_part(geom);
    }

    if (num_sides > 2)
      face_bond(&geom, &pyr);
    else {
      geom.del(FACES, vector<int>(1, 0));
      int apex = geom.add_vert(Vec3d(0, 0, ant.get_height(0) / 2 + ht));
      geom.add_face(0, 1, apex, -1);
      geom.add_face(1, 0, apex, -1);
    }
  }

  return Status::ok();
}

//--------------------------------------------------------------------

/*
void dipyramid::make_scal_part(Geometry &geom)
{
   double twist_ang = isnan(twist_angle) ? 0.0 : twist_angle;

   double edge = get_edge();
   double ring_ht = 0.0;
   if(subtype==subtype_scalenohedron)
      ring_ht = sin(twist_ang) * edge/2;

   double R1 = 0.0;
   double R2 = 0.0;
   if(subtype==subtype_scalenohedron)
      R1 = R2 = radius * cos(twist_ang);
   else if(subtype==subtype_dip_scalenohedron) {
      double ang = angle();
      R1 = edge*cos(ang/2 + twist_ang)/sin(ang);
      R2 = edge*cos(ang/2 - twist_ang)/sin(ang);
   }

   dipyramid dip(num_sides, step);
   dip.set_height(height);
   dip.set_twist_angle();
   dip.set_subtype(0);
   dip.make_poly_part(geom);
   vector<Vec3d> &verts = geom.raw_verts();
   for(int i=0; i<num_sides/2; i++) {
      verts[2*i]   = verts[2*i] * R1 + Vec3d(0,0, ring_ht);
      verts[2*i+1] = verts[2*i+1]*R2 + Vec3d(0,0,-ring_ht);
   }
}
*/

Status Polygon::make_dipyramid_scal_part(Geometry &geom)
{
  double rad = get_polygon_radius();
  double ht = (value_is_set(height[0])) ? height[0] : rad;
  double inrad = rad * cos(angle() / 2);
  double rad2 = value_is_set(radius[2]) ? radius[2] : inrad;

  Polygon pgon(num_sides, step);
  pgon.set_radius(0, rad);
  Geometry pg;
  pgon.add_polygon(pg);
  for (int i = 0; i < num_sides; i++) {
    geom.add_vert(pg.verts(i));
    if (num_sides > 2)
      geom.add_vert((rad2 / inrad) * 0.5 *
                    (pg.verts(i) + pg.verts((i + 1) % num_sides)));
    else
      geom.add_vert(Vec3d(0, (1 - 2 * i) * rad2, 0));
  }
  int apex_idx1 = geom.add_vert(Vec3d(0, 0, ht));
  int apex_idx2 = geom.add_vert(Vec3d(0, 0, -ht));

  for (int i = 0; i < 2 * num_sides; i++) {
    geom.add_face(i, (i + 1) % (2 * num_sides), apex_idx1, -1);
    geom.add_face((i + 1) % (2 * num_sides), i, apex_idx2, -1);
  }

  return Status::ok();
}

Status Polygon::make_dipyramid_part(Geometry &geom)
{
  double ht = (value_is_set(height[0])) ? height[0] : get_polygon_radius();
  if ((subtype == sub_default && value_is_set(twist_angle[0])) ||
      subtype == sub_dipyramid_trapezohedron) { // make trapezohedron
    if (num_sides == 2)
      return Status::error(
          "cannot make a trapezohedron from a digonal dipyramid");

    Polygon ant(*this);
    ant.set_type(antiprism);
    ant.set_subtype(sub_antiprism_trapezohedron);
    double twist_ang =
        value_is_set(twist_angle[0]) ? twist_angle[0] - angle() / 2 : 0.0;
    ant.set_twist_angle(0, twist_ang);
    double ant_ht = ht * (cos(twist_ang) / cos(angle() / 2) - 1);
    ant.set_height(0, ant_ht);
    ant.make_poly(geom);
    for (int i = 0; i < 2; i++) // make sure apex heights are correct
      geom.verts(geom.verts().size() - 1 - i) =
          Vec3d(0, 0, (1 - 2 * i) * (0.5 * ant_ht + ht));
  }
  else if (subtype == sub_dipyramid_scalenohedron)
    make_dipyramid_scal_part(geom);
  else { // make dipyramid types
    // Make the top part, possibly elongated or gyroelongated
    make_pyramid_part(geom);

    // Add the bottom pyramid
    Polygon pyr(*this);
    pyr.set_subtype(sub_default);
    pyr.set_twist_angle(0, NAN);
    Geometry pyr_geom;
    Status stat = pyr.make_pyramid_part(pyr_geom);
    if (stat.is_error())
      return stat;
    if (num_sides > 2)
      face_bond(&geom, &pyr_geom);
    else {
      geom.del(FACES, vector<int>(1, 0));
      int off = (subtype == sub_default) ? 0 : 2;
      int apex = geom.add_vert(Vec3d(0, 0, geom.verts(off)[2] - ht));
      geom.add_face(0 + off, 1 + off, apex, -1);
      geom.add_face(1 + off, 0 + off, apex, -1);
    }
  }

  return Status::ok();
}

//--------------------------------------------------------------------

static void prism_wrap(Geometry &geom, int num_wraps)
{
  vector<vector<int>> &faces = geom.raw_faces();
  int sz = faces.size();
  for (int i = 0; i < sz; i++) {
    int f_sz = faces[i].size();
    for (int w = 0; w < num_wraps - 1; w++) {
      if (i < 2) // double wind the caps
        faces[i].insert(faces[i].end(), faces[i].begin(),
                        faces[i].begin() + f_sz);
      else // repeat the side faces
        geom.add_face(faces[i]);
    }
  }
}

Status Polygon::make_cupola_part(Geometry &geom)
{
  if (subtype == sub_cupola_cuploid && !is_even(step * parts))
    return Status::error("cannot make a cupoloid unless the polygon "
                         "fraction, in lowest form, has an even denominator");

  double rad = get_polygon_radius();
  double ht;
  if (value_is_set(height[0]))
    ht = height[0];
  else if (value_is_set(edge[1])) {
    if (edge[1] - rad < -epsilon) // CHECK need fabs?
      return Status::error("edge too short to reach between polygons");

    double inrad1 = rad * cos(angle() / 2);
    double inrad2 =
        0.5 * get_polygon_edge() / tan(angle() / 4); // larger polygon
    double diff2 = pow(inrad2 - inrad1, 2);
    ht = (edge[1] * edge[1] > diff2) ? sqrt(edge[1] * edge[1] - diff2) : 0.0;
  }
  else
    ht = rad;

  bool even = is_even(step);
  Geometry cup_geom;
  vector<Vec3d> verts;
  vector<vector<int>> faces;
  int n2 = (2 - even) * num_sides;
  Polygon large(n2, step / (1 + even));
  large.set_edge(0, get_polygon_edge());
  large.add_polygon(cup_geom);
  if (even) {
    for (int i = 0; i < n2; i++)
      cup_geom.raw_faces()[0].push_back(i);
  }
  if (subtype == sub_cupola_cuploid)
    cup_geom.clear(FACES);

  cup_geom.transform(Trans3d::rotate(Vec3d::Z, -angle() / 4));

  int v_sz = cup_geom.verts().size();
  add_polygon(cup_geom, ht);
  int off = 0;
  for (int i = 0; i < 2 * num_sides; i++) {
    vector<int> face;
    int v = (i + n2 - off) % n2;
    face.push_back(v);
    face.push_back(v_sz + ((i + 1) / 2) % num_sides);
    if (is_even(i))
      face.push_back(v_sz + ((i / 2 + 1) % num_sides));
    face.push_back((v + 1) % n2);
    cup_geom.add_face(face);
  }

  if (subtype == sub_default || subtype == sub_cupola_cuploid) {
    geom.append(cup_geom);
  }
  else if (subtype == sub_cupola_elongated) {
    Polygon pri(large);
    pri.set_type(prism);
    pri.set_height(0, value_is_set(height[1]) ? height[1] : get_polygon_edge());
    pri.make_prism_part(geom);
    if (even)
      prism_wrap(geom, 2);
    face_bond(&geom, &cup_geom);
  }
  else if (subtype == sub_cupola_gyroelongated) {
    Polygon ant(large);
    if (value_is_set(height[1]))
      ant.set_height(0, height[1]);
    else
      ant.set_edge(1, get_polygon_edge());
    ant.make_antiprism_part(geom);
    if (even)
      prism_wrap(geom, 2);
    face_bond(&geom, &cup_geom);
  }

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_orthobicupola_part(Geometry &geom)
{
  // Make the top part, possibly elongated or gyroelongated
  Status stat = make_cupola_part(geom);
  if (stat.is_error())
    return stat;

  // Add the bottom cupola
  Polygon cup(*this);
  cup.set_subtype(sub_default);
  cup.set_twist_angle(0, NAN);
  Geometry cup_geom;
  cup.make_cupola_part(cup_geom);
  face_bond(&geom, &cup_geom, 0, 0, (subtype == sub_default)); // var'n in bases

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_gyrobicupola_part(Geometry &geom)
{
  // Make the top part, possibly elongated or gyroelongated
  Status stat = make_cupola_part(geom);
  if (stat.is_error())
    return stat;

  // Add the bottom cupola
  Polygon cup(*this);
  cup.set_subtype(sub_default);
  Geometry cup_geom;
  cup.make_cupola_part(cup_geom);
  face_bond(&geom, &cup_geom, 0, 0, (subtype != 0)); // variation in bases

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_snub_antiprism_part(Geometry &geom)
{
  const double sqrt_epsilon = sqrt(epsilon);
  const double ang_inc = step * M_PI / num_sides;
  const double s = sin(ang_inc);
  const double c = cos(ang_inc);
  double coeffs[5];
  coeffs[0] = (2 - 3 * c * c);
  coeffs[1] = 2 * (1 + c);
  coeffs[2] = 4 * c - 7;
  coeffs[3] = -4;
  coeffs[4] = 4;
  double sol[4];
  int num_roots = quartic(coeffs, sol);
  // fprintf(stderr, "\n%d radius values to test\n", num_roots);
  vector<double> valid;
  for (int i = 0; i < num_roots; i++) {
    double r = sol[i] / s;
    // fprintf(stderr, "\nradius (%d)\n   r = %.16f\n", i, r);
    double rt = 1 - pow(r - 0.5 / s, 2);
    // fprintf(stderr, "   h1^2 = %.16f\n", rt);
    if (rt < -sqrt_epsilon)
      continue;
    double h1 = rt > sqrt_epsilon ? sqrt(rt) : 0;
    rt = 3 / 4.0 - pow(r - 0.5 * c / s, 2);
    // fprintf(stderr, "   h2^2 = %.16f\n", rt);
    if (rt < -sqrt_epsilon)
      continue;
    double h2 = rt > sqrt_epsilon ? sqrt(rt) : 0;
    for (int j = 0; j < 2; j++) {
      h1 *= -1; // flip sign
      // fprintf(stderr, "      h1=%.16f\n", h1);
      // fprintf(stderr, "      h2=%.16f\n", h2);
      // fprintf(stderr, "      edge length = %.16f\n",
      // 2*r*r*(1-c)+(h1-h2)*(h1-h2));
      if (fabs(2 * r * r * (1 - c) + (h1 - h2) * (h1 - h2) - 1) > epsilon / s)
        continue;
      valid.push_back(r);
      valid.push_back(h1);
      valid.push_back(h2);
      if (fabs(h1) < epsilon || fabs(h2) < epsilon)
        break;
    }
  }

  int num_sols = valid.size() / 3;
  if (num_sols > 2)
    fprintf(stderr, "\n\n*** unexpected: snub-antiprism has more than 2 "
                    "solutions!!! (%d solutions)\n\n)",
            num_sols);

  if (subtype > num_sols - 1)
    return Status::error("subtype too large (the snub-antiprism only has %d "
                         "solution(s)",
                         num_sols);

  double r = valid[subtype * 3];
  double H =
      (valid[subtype * 3 + 1] + valid[3 * subtype + 2]) / 2; // new top layer ht
  double h = H - valid[subtype * 3 + 1]; // new upper inner layer height

  // fprintf(stderr, "\nFinal values\n   r = %.16f\n   H = %.16f\n   h =
  // %.16f\n",
  //      r, H, h);

  set_edge(0, 1.0);
  add_polygon(geom, H);
  vector<int> bond_face;
  for (int i = 0; i < num_sides; i++) {
    int a = num_sides + 2 * i;
    geom.add_vert(r * cos(2 * i * ang_inc) * Vec3d::X + h * Vec3d::Z +
                  -r * sin(2 * i * ang_inc) * Vec3d::Y);
    geom.add_vert(r * cos((2 * i + 1) * ang_inc) * Vec3d::X - h * Vec3d::Z +
                  -r * sin((2 * i + 1) * ang_inc) * Vec3d::Y);
    vector<int> face(3);
    face[0] = i;
    face[1] = a;
    face[2] = a + 1;
    geom.add_face(face);
    face[0] = i;
    face[1] = (i + 1) % num_sides;
    face[2] = a + 1;
    geom.add_face(face);
    face[2] = (i + 1) % num_sides;
    face[0] = a + 1;
    face[1] = num_sides + (2 * (i + 1)) % (2 * num_sides);
    geom.add_face(face);
    bond_face.push_back(a);
    bond_face.push_back(a + 1);
  }
  geom.add_face(bond_face);
  Geometry geom2 = geom;
  face_bond(&geom, &geom2, geom.faces().size() - 1, geom2.faces().size() - 1,
            1);
  // align dihedral axis with x-axis
  geom.transform(Trans3d::rotate(Vec3d::Z, M_PI * step / (2.0 * num_sides)));

  return Status::ok();
}

//--------------------------------------------------------------------

Status Polygon::make_crown_part(Geometry &geom)
{
  geom.clear_all();
  double twist_ang =
      value_is_set(get_twist_angle(0)) ? get_twist_angle(0) : 0.0;
  double e2;
  if (value_is_set(get_edge(1)))
    e2 = get_edge(1);
  else if (value_is_set(get_height(0))) {
    if (cos(twist_ang) == 0)
      return Status::error("height cannot be set for this twist angle");

    e2 = get_height(0) / cos(twist_ang);
  }
  else
    e2 = get_polygon_radius();

  Geometry pgon[2];
  add_polygon(pgon[0]);
  add_polygon(pgon[1]);
  pgon[1].transform(Trans3d::rotate(0, 0, -angle() / 2));
  int num = 4 * num_sides;
  for (int i = 0; i < 2 * num_sides; i++) {
    int p_idx = i % 2;
    int v_idx = i / 2;
    Vec3d v =
        Trans3d::rotate(pgon[p_idx].verts(v_idx), (1 - 2 * p_idx) * twist_ang) *
        Vec3d(0, 0, e2 / 2);
    geom.add_vert(pgon[p_idx].verts(v_idx) + v);
    geom.add_vert(pgon[p_idx].verts(v_idx) - v);
    geom.add_face(2 * i, 2 * i + 1, (2 * (i + 1)) % num,
                  (2 * (i + 1) + 1) % num, -1);
    geom.add_face((2 * (i + 0) + p_idx) % num, (2 * (i + 1) + !p_idx) % num,
                  (2 * (i + 3) + p_idx) % num, (2 * (i + 2) + !p_idx) % num,
                  -1);
  }

  return Status::ok();
}

Status Polygon::make_crown_full(Geometry &geom)
{
  int N = num_sides * parts;
  int step1 = (((step * parts) % N) + N) % N;
  int step2 = !value_is_set(get_twist_angle(1))
                  ? 1
                  : int(floor(rad2deg(get_twist_angle(1)) + 0.5));
  bool is_prism = is_even(step1 + step2);

  if (is_prism) {
    Polygon pri(N);
    pri.set_radius(0, get_polygon_radius());
    double height = value_is_set(get_height(0)) ? get_height(0) : get_edge(1);
    pri.set_height(0, height);
    Status stat = pri.make_prism_part(geom);
    if (stat.is_error())
      return stat;
    geom.clear(FACES);
    int offset = (step2 - step1) / 2;
    for (int i = 0; i < N; i++) {
      geom.add_face((i - offset + N) % N, (i + step1) % N + N,
                    (i + step1 + offset + N) % N, (i) % N + N, -1);

      geom.add_face((i - offset + N) % N + N, (i + step1) % N,
                    (i + step1 + offset + N) % N + N, (i) % N, -1);
    }
  }
  else {
    Polygon ant(N);
    ant.set_radius(0, get_polygon_radius());
    ant.set_height(0, get_height(0));
    ant.set_edge(1, get_edge(1));
    Status stat = ant.make_antiprism_part(geom);
    if (stat.is_error())
      return stat;
    geom.clear(FACES);
    bool odd1 = !is_even(step1);
    bool odd2 = !is_even(step2);
    for (int i = 0; i < N; i++) {
      geom.add_face(
          (i) % N,
          (i - 1 + (step2 + odd2) / 2 + (step1 + odd1) / 2 + N) % N + N,
          (i + step2) % N,
          (i - 1 + (step2 + odd2) / 2 - (step1 - odd1) / 2 + N) % N + N, -1);
      geom.add_face(
          (i) % N,
          (i - 1 + (step1 + odd1) / 2 + (step2 + odd2) / 2 + N) % N + N,
          (i + step1) % N,
          (i - 1 + (step1 + odd1) / 2 - (step2 - odd2) / 2 + N) % N + N, -1);
    }
  }

  return Status::ok();
}

} // namespace anti
