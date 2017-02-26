/*
   Copyright (c) 2010-2016, Roger Kaufman

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

/*!\file normals.cc
 * \brief classes for working with normals
*/

#include <string>
#include <vector>

#include "geometryutils.h"

using std::string;
using std::vector;

namespace anti {

// face normal
Normal::Normal(const Geometry &geom, const int face_idx, Vec3d C, double eps)
{
  const vector<int> &face = geom.faces()[face_idx];
  const vector<Vec3d> &verts = geom.verts();

  normal = face_norm(verts, face);

  if (!C.is_set())
    C = centroid(verts);

  double D = vdot(verts[face[0]] - C, normal);

  // test case centroid is on the hemi's face
  if (double_eq(vdot(centroid(verts, face) - C, normal), 0.0, eps))
    direction = 0;
  else
    direction = double_compare(D, 0.0, eps);
}

// edge and vertex normals which have already been calculated from FaceNormals
Normal::Normal(const Geometry &geom, const Vec3d &norm, const int v_idx,
               Vec3d C, double eps)
{
  const vector<Vec3d> &verts = geom.verts();

  normal = norm;

  if (!C.is_set())
    C = centroid(verts);

  double D = vdot(verts[v_idx] - C, normal);

  direction = double_compare(D, 0.0, eps);
}

void FaceNormals::refresh(const Geometry &geom, Vec3d cent, double eps)
{
  ngeom = &geom;

  if (!cent.is_set())
    cent = centroid((*ngeom).verts());

  normals.clear();
  for (unsigned int i = 0; i < geom.faces().size(); i++)
    normals.push_back(Normal(*ngeom, i, cent, eps));
}

// private common code for averaging edge and vertex normals
Vec3d FaceNormals::average_normals(vector<int> &face_idx,
                                   const string &average_pattern)
{
  Vec3d norm;

  if (face_idx.size()) {
    vector<Vec3d> f_normals;
    for (int i : face_idx) {
      Vec3d normal = normals[i].raw(); // default
      for (char j : average_pattern) {
        if (j == 'r')
          normal = normals[i].raw();
        else if (j == 'o')
          normal = normals[i].outward();
        else if (j == 'i')
          normal = normals[i].inward();
        else if (j == 'u')
          normal.to_unit();
      }
      f_normals.push_back(normal);
    }

    norm = centroid(f_normals);
  }

  return norm;
}

// the edge normal is centroid of all the face normals of which edge
// is a part of those faces
Normal FaceNormals::edge_normal(int v_idx1, int v_idx2,
                                const string &average_pattern, Vec3d cent,
                                double eps)
{
  vector<int> edge = make_edge(v_idx1, v_idx2);

  vector<int> face_idx = find_faces_with_edge((*ngeom).faces(), edge);

  Vec3d norm = average_normals(face_idx, average_pattern);

  return (norm.is_set() ? Normal(*ngeom, norm, v_idx1, cent, eps) : Normal());
}

// the vertex normal is centroid of all the face normals of which vertex
// is a part of those faces
Normal FaceNormals::vert_normal(int v_idx, const string &average_pattern,
                                Vec3d cent, double eps)
{
  vector<int> face_idx = find_faces_with_vertex((*ngeom).faces(), v_idx);

  Vec3d norm = average_normals(face_idx, average_pattern);

  return (norm.is_set() ? Normal(*ngeom, norm, v_idx, cent, eps) : Normal());
}

} // namespace anti
