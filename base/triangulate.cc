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

/* \file triangulate.cc
   convert polyhedron faces to triangles
*/

#include <algorithm>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "geometry.h"
#include "tesselator/glu.h"

using std::map;
using std::vector;

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#define APIENTRY

namespace anti {

void triangulate_basic(Geometry &geom, bool sq_diag, Color inv,
                       vector<int> *fmap)
{
  vector<int> del_faces;
  int idx;
  unsigned int num_faces = geom.faces().size();
  if (fmap)
    fmap->clear();
  for (unsigned int i = 0; i < num_faces; ++i) {
    // Store the index of the first triangle part of this face
    if (fmap)
      fmap->push_back(geom.faces().size() - num_faces);
    int fsz = geom.faces(i).size();
    if (fsz <= 3)
      continue;
    Color col;
    col = geom.colors(FACES).get((int)i);
    del_faces.push_back(i);
    vector<int> tface(3);
    if (fsz == 4 && sq_diag) {
      bool st = (geom.face_v(i, 0) - geom.face_v(i, 2)).len2() >
                (geom.face_v(i, 1) - geom.face_v(i, 3)).len2();
      tface[0] = geom.faces_mod(i, 3 + st);
      tface[1] = geom.faces(i, 0 + st);
      tface[2] = geom.faces(i, 1 + st);
      idx = geom.add_face(tface);
      geom.colors(FACES).set(idx, col);
      tface[0] = geom.faces(i, 1 + st);
      tface[1] = geom.faces(i, 2 + st);
      tface[2] = geom.faces_mod(i, 3 + st);
      idx = geom.add_face(tface);
      geom.colors(FACES).set(idx, col);
      if (inv.is_set())
        geom.colors(EDGES).set(geom.add_edge(tface[0], tface[2]), inv);
      continue;
    }

    tface[2] = geom.add_vert(geom.face_cent(i));
    geom.colors(VERTS).set(tface[2], inv);
    for (int j = 0; j < fsz; ++j) {
      tface[0] = geom.faces(i, j);
      tface[1] = geom.faces_mod(i, j + 1);
      idx = geom.add_face(tface);
      geom.colors(FACES).set(idx, col);
      if (inv.is_set())
        geom.colors(EDGES).set(geom.add_edge(tface[0], tface[2]), inv);
    }
  }
  geom.del(FACES, del_faces);
}

struct face_tris {
  Geometry *geom;
  Color col;
  Color inv;
  vector<int> idxs;
  int e_idx;
  vector<Vec3d *> extra_verts;
  vector<int *> extra_idxs;

  face_tris(Geometry *geo, Color c, Color i)
      : geom(geo), col(c), inv(i), e_idx(-1)
  {
  }
  ~face_tris()
  {
    update_faces();
    for (auto &extra_vert : extra_verts)
      delete extra_vert;
    for (auto &extra_idx : extra_idxs)
      delete extra_idx;
  }
  void update_faces();
};

void face_tris::update_faces()
{
  map<vector<int>, int> edge_cnts;
  for (unsigned int i = 0; i < idxs.size() / 3; ++i) {
    vector<int> face(3);
    for (int j = 0; j < 3; ++j)
      face[j] = idxs[i * 3 + j];
    int idx = geom->add_face(face);
    geom->colors(FACES).set(idx, col);

    if (inv.is_set())
      for (int j = 0; j < 3; ++j)
        edge_cnts[make_edge(face[j], face[(j + 1) % 3])]++;
  }

  if (inv.is_set()) {
    const vector<vector<int>> &edges = geom->edges();
    map<vector<int>, int>::const_iterator mi;
    for (mi = edge_cnts.begin(); mi != edge_cnts.end(); ++mi) {
      if (mi->second == 2) { // new edge internal to a face
        if (find(edges.begin(), edges.end(), mi->first) == edges.end()) {
          // new edge is not an explicit edge so add as an invisible edge
          int eidx = geom->add_edge(mi->first);
          geom->colors(EDGES).set(eidx, inv);
        }
      }
    }
  }
}

// Dummy callback ensures localGL_TRIANGLES are used
extern "C" APIENTRY void tri_eflag(localGLboolean) {}

extern "C" APIENTRY void tri_vert(void *vdata, void *data)
{
  int idx = *(int *)vdata;
  face_tris *f_tris = (face_tris *)data;
  f_tris->idxs.push_back(idx);
}

extern "C" APIENTRY void tri_combine(localGLdouble coords[3], localGLdouble **,
                                     localGLfloat *, void **dataOut, void *data)
{
  auto *v = new Vec3d(coords[0], coords[1], coords[2]);
  face_tris *f_tris = (face_tris *)data;
  auto *pidx = new int(f_tris->geom->add_vert(*v));
  f_tris->geom->colors(VERTS).set(*pidx, f_tris->inv);
  *dataOut = pidx;
  f_tris->extra_verts.push_back(v);
  f_tris->extra_idxs.push_back(pidx);
}

class anti_tesselator {
private:
  localGLUtesselator *tess;

public:
  anti_tesselator();
  ~anti_tesselator();
  bool set_winding_rule(unsigned int winding_rule);
  operator localGLUtesselator *() { return tess; }
};

anti_tesselator::anti_tesselator()
{
  tess = localgluNewTess();
  set_winding_rule(TESS_WINDING_NONZERO);
  localgluTessCallback(tess, localGLU_TESS_EDGE_FLAG,
                       (_localGLfuncptr)tri_eflag);
  localgluTessCallback(tess, localGLU_TESS_VERTEX_DATA,
                       (_localGLfuncptr)tri_vert);
  localgluTessCallback(tess, localGLU_TESS_COMBINE_DATA,
                       (_localGLfuncptr)tri_combine);
}

anti_tesselator::~anti_tesselator() { localgluDeleteTess(tess); }

bool anti_tesselator::set_winding_rule(unsigned int winding_rule)
{
  if (winding_rule < TESS_WINDING_ODD ||
      winding_rule > TESS_WINDING_ABS_GEQ_TWO)
    return false;

  localgluTessProperty(tess, localGLU_TESS_WINDING_RULE, winding_rule);
  return true;
}

void triangulate(Geometry &geom, Color inv, unsigned int winding,
                 vector<int> *fmap)
{
  anti_tesselator tess;
  tess.set_winding_rule(winding);
  vector<vector<int>> faces = geom.faces();
  vector<vector<int>> impl_edges;
  geom.get_impl_edges(impl_edges);
  map<int, Color> fcolmap;
  fcolmap = geom.colors(FACES).get_properties();
  geom.clear(FACES);

  const vector<Vec3d> &verts = geom.verts();
  if (fmap)
    fmap->clear();
  for (unsigned int i = 0; i < faces.size(); i++) {
    // Store the index of the first triangle part of this face
    if (fmap)
      fmap->push_back(geom.faces().size());
    if (faces[i].size() < 3)
      continue;
    Color col;
    auto mi = fcolmap.find(i);
    if (mi != fcolmap.end())
      col = mi->second;

    face_tris f_tris(&geom, col, inv);
    localgluTessBeginPolygon(tess, &f_tris);
    for (int &j : faces[i]) {
      double vtx[3]; // tesselator sometimes fails when using doubles (?)
      vtx[0] = (float)verts[j][0];
      vtx[1] = (float)verts[j][1];
      vtx[2] = (float)verts[j][2];
      localgluTessVertex(tess, vtx, &j);
    }
    localgluTessEndPolygon(tess);
  }
}

} // namespace anti
