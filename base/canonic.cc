/*
   Copyright (c) 2003-2016, Adrian Rossiter, Roger Kaufman
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
   Name: canonic.cc
   Description: canonicalize a polyhedron
                Implementation of George Hart's canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
   Project: Antiprism - http://www.antiprism.com
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "geometry.h"
#include "geometryinfo.h"

using std::string;
using std::vector;

namespace anti {

// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
void canonicalize_mm(Geometry *geom, double edge_factor, double plane_factor,
                     int n, int divergence_test, int rep_count,
                     bool planar_only, double eps)
{
  // do a scale to get edges close to 1
  GeometryInfo info(*geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    geom->transform(Trans3d::scale(1 / scale));

  const vector<Vec3d> &verts = geom->verts();
  const vector<vector<int>> &faces = geom->faces();
  vector<vector<int>> edges;
  geom->get_impl_edges(edges);

  GeometryInfo rep(*geom);
  rep.set_center(geom->centroid());
  double starting_radius = rep.vert_dist_lims().max;

  double max_diff2 = 0;
  Vec3d origin(0, 0, 0);
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)n;) {
    vector<Vec3d> old_verts = verts;

    if (!planar_only) {
      vector<Vec3d> near_pts;
      for (auto &edge : edges) {
        Vec3d P = geom->edge_nearpt(edge, origin);
        near_pts.push_back(P);
        Vec3d offset = edge_factor * (P.len() - 1) * P;
        geom->verts(edge[0]) -= offset;
        geom->verts(edge[1]) -= offset;
      }

      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        geom->verts(i) -= cent_near_pts;
    }

    // Make a copy of verts. zero out.
    vector<Vec3d> vs = verts;
    for (auto &v : vs)
      v = Vec3d(0, 0, 0);

    for (unsigned int ff = cnt; ff < faces.size() + cnt; ff++) {
      int f = ff % faces.size();
      if (faces[f].size() == 3)
        continue;
      Vec3d norm = geom->face_norm(f).unit();
      Vec3d f_cent = geom->face_cent(f);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;
      for (int v : faces[f])
        vs[v] += vdot(plane_factor * norm, f_cent - verts[v]) * norm;
    }

    // adjust vertices post-loop
    for (unsigned int i = 0; i < vs.size(); i++)
      geom->verts(i) += vs[i];

    max_diff2 = 0;
    for (unsigned int i = 0; i < verts.size(); i++) {
      double diff2 = (verts[i] - old_verts[i]).len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // increment count here for reporting
    cnt++;

    if (rep_count > 0 && (cnt) % rep_count == 0)
      fprintf(stderr, "%-15d max_diff=%12.10g\n", cnt, sqrt(max_diff2));

    if (sqrt(max_diff2) < eps)
      break;

    // see if radius is expanding or contracting unreasonably
    if (divergence_test > 0) {
      rep.set_center(geom->centroid());
      if ((rep.vert_dist_lims().max > starting_radius * divergence_test) ||
          (rep.vert_dist_lims().max < starting_radius / divergence_test)) {
        fprintf(stderr, "Probably Diverging. Breaking out.\n");
        break;
      }
    }
  }
  if (rep_count > -1) {
    fprintf(stderr, "\n%-15d final max_diff=%12.10g\n", cnt, sqrt(max_diff2));
    fprintf(stderr, "\n");
  }
}

vector<Vec3d> reciprocalN(Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<Vec3d> normals;
  for (const auto &face : faces) {
    Vec3d centroid(0, 0, 0);
    Vec3d normal(0, 0, 0);
    double avgEdgeDist = 0;

    int v1 = face.at(face.size() - 2);
    int v2 = face.at(face.size() - 1);
    for (int v3 : face) {
      centroid += verts[v3];
      // normal += orthogonal(verts[v1], verts[v2], verts[v3]);
      normal += vcross(verts[v3] - verts[v2],
                       verts[v2] - verts[v1]); // orthogonal without call
      // avgEdgeDist += tangentPoint(verts[v1], verts[v2]).len();
      Vec3d d = verts[v2] - verts[v1];
      // prevent division by zero
      // avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.len2()) * d)).len();
      // // tangentPoint without call
      double vdt;
      if (d[0] == 0 && d[1] == 0 && d[2] == 0)
        vdt = 0;
      else
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len(); // tangentPoint without call
      v1 = v2;
      v2 = v3;
    }
    centroid *= 1.0 / face.size();
    normal.to_unit();
    avgEdgeDist /= face.size();

    // reciprocal call replace below:
    // Vec3d ans = reciprocal(normal * vdot(centroid,normal));
    Vec3d v = normal * vdot(centroid, normal);
    // prevent division by zero
    Vec3d ans;
    if (v[0] == 0 && v[1] == 0 && v[2] == 0)
      ans = v;
    else {
      ans = v * 1.0 / v.len2();
      ans *= (1 + avgEdgeDist) / 2;
    }
    normals.push_back(ans);
  }

  return normals;
}

vector<Vec3d> reciprocalC(Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len2();
  return centers;
}

vector<Vec3d> reciprocalC_len(Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len();
  return centers;
}

// Addition to algorithm by Adrian Rossiter. Finds the correct centroid for the
// canonical
Vec3d edge_nearpoints_centroid(Geometry &geom, Vec3d cent)
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  int e_sz = edges.size();
  Vec3d e_cent(0, 0, 0);
  for (int e = 0; e < e_sz; ++e)
    e_cent += geom.edge_nearpt(edges[e], cent);
  return e_cent / double(e_sz);
}

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
void canonicalize_cn(Geometry *geom, int n, char method, int divergence_test,
                     int rep_count, double eps)
{
  // do a scale to get edges close to 1
  GeometryInfo info(*geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    geom->transform(Trans3d::scale(1 / scale));

  Geometry dual;
  get_dual(&dual, *geom, 0);
  dual.clear_cols();
  const vector<Vec3d> &verts = geom->verts();
  const vector<Vec3d> &d_verts = dual.verts();

  GeometryInfo rep(*geom);
  rep.set_center(geom->centroid());
  double starting_radius = rep.vert_dist_lims().max;

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)n;) {
    vector<Vec3d> old_verts = verts;
    Vec3d e_cent;

    switch (method) {
    // base/dual canonicalize method
    case 'n':
      dual.raw_verts() = reciprocalN(*geom);
      geom->raw_verts() = reciprocalN(dual);
      e_cent = edge_nearpoints_centroid(*geom, Vec3d(0, 0, 0));
      geom->transform(Trans3d::transl(-0.1 * e_cent));
      break;

    // adjust vertices with side effect of planarization. len2() version
    case 'p':
      // move centroid to origin for balance
      dual.raw_verts() = reciprocalC(*geom);
      geom->transform(Trans3d::transl(-centroid(d_verts)));
      geom->raw_verts() = reciprocalC(dual);
      geom->transform(Trans3d::transl(-centroid(verts)));
      break;

    // adjust vertices with side effect of planarization. len() version
    case 'q':
      // move centroid to origin for balance
      dual.raw_verts() = reciprocalC_len(*geom);
      geom->transform(Trans3d::transl(-centroid(d_verts)));
      geom->raw_verts() = reciprocalC_len(dual);
      geom->transform(Trans3d::transl(-centroid(verts)));
      break;

    case 'x':
      geom->face_cents(dual.raw_verts());
      dual.face_cents(geom->raw_verts());
      break;
    }

    max_diff2 = 0;
    for (unsigned int i = 0; i < verts.size(); i++) {
      double diff2 = (verts[i] - old_verts[i]).len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // increment count here for reporting
    cnt++;

    if (rep_count > 0 && (cnt) % rep_count == 0)
      fprintf(stderr, "%-15d max_diff=%12.10g\n", cnt, sqrt(max_diff2));

    if (sqrt(max_diff2) < eps)
      break;

    // see if radius is expanding or contracting unreasonably
    if (divergence_test > 0) {
      rep.set_center(geom->centroid());
      if ((rep.vert_dist_lims().max > starting_radius * divergence_test) ||
          (rep.vert_dist_lims().min < starting_radius / divergence_test)) {
        fprintf(stderr, "Probably Diverging. Breaking out.\n");
        break;
      }
    }
  }

  if (rep_count > -1) {
    fprintf(stderr, "\n%-15d final max_diff=%12.10g\n", cnt, sqrt(max_diff2));
    fprintf(stderr, "\n");
  }
}

} // namespace anti
