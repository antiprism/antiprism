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

/* \file c_hull.cc
   \brief wrapper around qhull for convex hulls and delaunay edges.
*/

#include "boundbox.h"
#include "geometry.h"
#include "geometryutils.h"
#include "mathutils.h"
#include "qhull/qhull_ra.h"
#include "utils.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <string>
#include <vector>

#define NUM_TO_USE_DELAUNAY 20

using std::map;
using std::pair;
using std::string;
using std::swap;
using std::vector;

namespace anti {

static void qhull_cleanup(qhT *qh)
{
  qh_freeqhull(qh, !qh_ALL); // free long memory
  int curlong, totlong;
  qh_memfreeshort(qh, &curlong, &totlong); // free short mem and mem allocator
}

static bool make_hull(Geometry &geom, bool append, string qh_args, char *errmsg)
{
  vector<Vec3d> verts = geom.verts();
  Vec3d cent = geom.centroid();

  map<int, Color> vcols;
  vcols = geom.colors(VERTS).get_properties();

  const int dim = 3;
  auto *points = new coordT[verts.size() * dim];

  for (unsigned i = 0; i < verts.size(); i++)
    for (int j = 0; j < dim; j++)
      points[i * dim + j] = verts[i][j];

  qh_args.insert(0, "qhull o ");

  boolT ismalloc = False;  // don't free points in qh_freeqhull() or realloc
  FILE *outfile = nullptr; // suppress output from qh_produce_output()
  FILE *errfile = fopen("/dev/null", "w"); // suppress qhull error messages
  if (!errfile)                            // must be a valid pointer
    errfile = stderr;

  qhT qh_val;
  qhT *qh = &qh_val;
  QHULL_LIB_CHECK
  qh_zero(qh, errfile);

  int ret = qh_new_qhull(qh, dim, verts.size(), points, ismalloc,
                         (char *)qh_args.c_str(), outfile, errfile);

  if (ret) {
    if (errmsg)
      snprintf(errmsg, MSG_SZ, "error calculating convex hull");
    qhull_cleanup(qh);
    delete[] points;
    return false;
  }

  map<size_t, int> vert_order;
  vertexT *vertex;

  if (!append) {
    geom.clear_all();

    int i = 0;
    FORALLvertices
    {
      size_t idx = (vertex->point - points) / dim;
      vert_order[idx] = i++;
      int v_idx = geom.add_vert(verts[idx]);
      geom.colors(VERTS).set(v_idx, vcols[idx]);
    }
  }

  facetT *facet;
  FORALLfacets
  {
    vector<int> face, ordered_face;
    vertexT *vid; //, **vidp;
    int vid_i, vid_n;
    FOREACHsetelement_i_(qh, vertexT, facet->vertices, vid)
    {
      if (!append)
        face.push_back(vert_order[(vid->point - points) / dim]);
      else
        face.push_back((vid->point - points) / dim);
    }
    if (face.size() > 3) {
      // order the face vertices so joining them sequentially will
      // form the polygon
      unsigned int r_cnt = 0;
      vector<int> lns;
      ridgeT *ridge; //, **ridgep;
      int ridge_i, ridge_n;
      FOREACHsetelement_i_(qh, ridgeT, facet->ridges, ridge)
      {
        r_cnt++;
        FOREACHsetelement_i_(qh, vertexT, ridge->vertices, vid)
        {
          if (!append)
            lns.push_back(vert_order[(vid->point - points) / dim]);
          else
            lns.push_back((vid->point - points) / dim);
        }
      }
      if (lns.size() != 2 * r_cnt) {
        snprintf(errmsg, MSG_SZ, "ridge included more than two points");
      }

      ordered_face.push_back(lns[0]);
      int pt = lns[1];
      for (unsigned int j = 1; j < r_cnt; j++) {
        ordered_face.push_back(pt);
        for (unsigned int k = j; k < r_cnt; k++) {
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
    else {
      ordered_face = face;
    }
    int f_no = geom.add_face(ordered_face);
    if (vdot(geom.face_norm(f_no), geom.face_v(f_no, 0) - cent) < -epsilon)
      reverse(geom.raw_faces()[f_no].begin(), geom.raw_faces()[f_no].end());
  }

  qhull_cleanup(qh);
  delete[] points;

  return true;
}

static int dimension_safe_make_hull(Geometry &geom, bool append, string qh_args,
                                    char *errmsg)
{
  // an empty geom should be the only reason an error can occur
  if (!(geom.verts().size())) {
    if (errmsg)
      snprintf(errmsg, MSG_SZ, "convex hull could not be created. no vertices");
    return -1;
  }

  int dimension = 3;
  if (!make_hull(geom, append, qh_args, errmsg)) {
    // if make_hull fails
    // assume point, line, or polygon 3 to limit

    BoundBox bb(geom.verts());
    Vec3d min = bb.get_min();
    Vec3d max = bb.get_max();
    double D = bb.max_width();

    // store original edge number. store edge widths. check to see if widths are
    // zero
    vector<pair<double, int>> e(3);
    for (unsigned int i = 0; i < 3; i++) {
      e[i].second = i;
      e[i].first = max[i] - min[i];
    }

    // sort on width
    sort(e.begin(), e.end());

    // check for special cases for line or point. if so, then no make_hull()
    // needed
    if (!e[0].first && !e[1].first) {
      dimension = 1;
      if (!e[2].first)
        dimension = 0;
    }

    // if dimension is still set at 3, no special case was found so try finding
    // polygon or line
    if (dimension == 3) {
      dimension = 2;

      // create the parallel vectors
      // shortest width (could be zero width) is the first parallel
      // next narrowest width is the second parallel
      vector<Vec3d> pvec(2);
      for (unsigned int i = 0; i < 2; i++) {
        pvec[i] = Vec3d(0, 0, 0);
        pvec[i][e[i].second] = D;
      }

      // make first point in a direction from the centre that is parallel
      // to the first edge (also narrowest edge) in the sort list
      // add the point and keep track of it
      Vec3d cent = bb.get_centre();
      Vec3d point1 = cent + pvec[0];
      geom.add_vert(point1);

      // 2D test for polygon
      if (!make_hull(geom, append, qh_args, errmsg)) {
        // if test fails need to test with another point to verify not 2D
        // remove point 1 and enter point 2 (also second narrowest edge) in the
        // sort list
        // add the point and keep track of it

        // delete point1
        int v_idx = find_vert_by_coords(geom, point1, DBL_MIN);
        if (v_idx != -1)
          geom.del(VERTS, v_idx);

        Vec3d point2 = cent + pvec[1];
        geom.add_vert(point2);

        if (!make_hull(geom, append, qh_args, errmsg)) {
          // if still error then vertices are on a line. Add a second point
          dimension = 1;

          // re-enter point 1 so there are 2 extra points to verify a line
          // If fails then 1D
          geom.add_vert(point1);

          // 1D test for line
          if (!make_hull(geom, append, qh_args, errmsg)) {
            // 0 dimensional geom should not have gotten in here. instead, make
            // this an error so it can be spotted
            if (errmsg)
              snprintf(errmsg, MSG_SZ,
                       "convex hull failed even after checking for a line");
            return -1;
          }

          // delete point1
          v_idx = find_vert_by_coords(geom, point1, DBL_MIN);
          if (v_idx != -1)
            geom.del(VERTS, v_idx);
        }

        // delete point2
        v_idx = find_vert_by_coords(geom, point2, DBL_MIN);
        if (v_idx != -1)
          geom.del(VERTS, v_idx);
      }

      // delete point1
      // if dimension is 1, point1 has already been deleted
      if (dimension != 1) {
        int v_idx = find_vert_by_coords(geom, point1, DBL_MIN);
        if (v_idx != -1)
          geom.del(VERTS, v_idx);
      }
    }

    // for dimension less than 2, what follows may have never gone through
    // make_hull()
    // for points and line strip any pre-existing faces or edges (that probably
    // won't exist)
    if (dimension < 2 && !append) {
      geom.clear(FACES);
      geom.clear(EDGES);
    }

    // if dimension is 0, if append is false, reduce to one point
    if (dimension == 0) {
      if (!append)
        merge_coincident_elements(geom, "v", epsilon);
    }
    // what should be left will be the line or polygon
    // if 1 dimensional, add one edge between the two points
    // a point may be double at the end so the vertices have to be sorted so
    // first and last are the end points for the edge connection
    else if (dimension == 1) {
      merge_coincident_elements(geom, "s", epsilon);
      // if not append, make sure all that exist is the first and last point
      if (!append) {
        vector<Vec3d> verts = geom.verts();
        Vec3d vert1 = verts[0];
        Vec3d vert2 = verts[geom.verts().size() - 1];
        geom.clear_all();
        geom.add_vert(vert1);
        geom.add_vert(vert2);
      }
      geom.add_edge(0, geom.verts().size() - 1);
    }
  }

  return dimension;
}

Status add_hull(Geometry &geom, string qh_args, int *dim)
{
  Status stat;
  char errmsg[MSG_SZ];
  int ret = dimension_safe_make_hull(geom, true, qh_args, errmsg);
  if (dim)
    *dim = ret;
  if (ret < 0)
    stat.set_error(errmsg);
  return stat;
}

Status set_hull(Geometry &geom, string qh_args, int *dim)
{
  Status stat;
  char errmsg[MSG_SZ];
  int ret = dimension_safe_make_hull(geom, false, qh_args, errmsg);
  if (dim)
    *dim = ret;
  if (ret < 0) {
    stat.set_error(errmsg);
    geom.clear_all();
  }
  return stat;
}

Status get_delaunay_edges(const vector<Vec3d> &verts,
                          map<pair<int, int>, int> &edges, string qh_args)
{
  const int dim = 3;
  auto *points = new coordT[verts.size() * dim];

  for (unsigned i = 0; i < verts.size(); i++)
    for (int j = 0; j < dim; j++)
      points[i * dim + j] = verts[i][j];

  qh_args.insert(0, "qhull d Qbb QJ o ");

  boolT ismalloc = False;  // don't free points in qh_freeqhull() or realloc
  FILE *outfile = nullptr; // suppress output from qh_produce_output()
  FILE *errfile = fopen("/dev/null", "w"); // suppress qhull error messages
  if (!errfile)                            // must be a valid pointer
    errfile = stderr;

  qhT qh_val;
  qhT *qh = &qh_val;
  QHULL_LIB_CHECK
  qh_zero(qh, errfile);

  if (qh_new_qhull(qh, dim, verts.size(), points, ismalloc,
                   (char *)qh_args.c_str(), outfile, errfile)) {
    qhull_cleanup(qh);
    delete[] points;
    return Status::error("error calculating delaunay triangulation");
  }

  edges.clear();
  facetT *facet;
  FORALLfacets
  {
    if (!facet->good || facet->upperdelaunay)
      continue;
    vector<int> tet;
    vertexT *vid;
    int vid_i, vid_n;
    FOREACHsetelement_i_(qh, vertexT, facet->vertices, vid)
    {
      tet.push_back((vid->point - qh->first_point) / (dim + 1));
    }

    if (tet.size() < 4)
      continue;

    bool coplanar = true;
    Vec3d v1 = verts[tet[1]] - verts[tet[0]];
    Vec3d v2 = verts[tet[2]] - verts[tet[0]];
    for (unsigned int i = 3; i < tet.size(); i++)
      if (fabs(vtriple(v1, v2, verts[tet[i]] - verts[tet[0]])) > epsilon) {
        coplanar = false;
        break;
      }
    if (coplanar)
      continue;

    sort(tet.begin(), tet.end());
    for (unsigned int i = 0; i < tet.size() - 1; i++)
      for (unsigned int j = i + 1; j < tet.size(); j++) {
        edges[pair<int, int>(tet[i], tet[j])]++;
      }
  }

  qhull_cleanup(qh);
  delete[] points;

  return Status::ok();
}

double get_min_vert_to_vert_dist(const vector<Vec3d> &verts, double sig_dist)
{
  double sig_dist2 = sig_dist * sig_dist;
  double min_sig_dist = 1e100;
  double mag2;

  if (!verts.size())
    min_sig_dist = 0;
  else if (verts.size() < NUM_TO_USE_DELAUNAY) {
    Vec3d dvec;
    for (unsigned int i = 0; i < verts.size() - 1; i++)
      for (unsigned int j = i + 1; j < verts.size(); j++) {
        mag2 = (verts[j] - verts[i]).len2();
        if (mag2 >= sig_dist2 && mag2 < min_sig_dist)
          min_sig_dist = mag2;
      }
  }
  else {
    map<pair<int, int>, int> edges;
    if (!get_delaunay_edges(verts, edges, ""))
      return -1; // error

    map<pair<int, int>, int>::iterator mi;
    for (mi = edges.begin(); mi != edges.end(); ++mi) {
      mag2 = (verts[mi->first.second] - verts[mi->first.first]).len2();
      if (mag2 >= sig_dist2 && mag2 < min_sig_dist)
        min_sig_dist = mag2;
    }
  }

  return (min_sig_dist < 1e99) ? sqrt(min_sig_dist) : sig_dist;
}

Status get_voronoi_cells(const vector<Vec3d> &verts, vector<Geometry> *cells,
                         string qh_args)
{
  const int dim = 3;
  auto *points = new coordT[verts.size() * dim];

  for (unsigned i = 0; i < verts.size(); i++)
    for (int j = 0; j < dim; j++)
      points[i * dim + j] = verts[i][j];

  qh_args.insert(0, "qhull v o ");

  boolT ismalloc = False;  /* don't free points in qh_freeqhull() or realloc*/
  FILE *outfile = nullptr; /* output from qh_produce_output()   nullptr ?*/
  FILE *errfile = fopen("/dev/null", "w"); // suppress qhull error messages
  // FILE *errfile= nullptr;    /* error messages from qhull code */

  qhT qh_val;
  qhT *qh = &qh_val;
  QHULL_LIB_CHECK
  qh_zero(qh, errfile);

  if (qh_new_qhull(qh, dim, verts.size(), points, ismalloc,
                   (char *)qh_args.c_str(), outfile, errfile)) {
    qhull_cleanup(qh);
    delete[] points;
    return Status::error("error calculating voronoi cells");
  }

  // qh_setvoronoi_all();
  qh_clearcenters(qh, qh_ASvoronoi);
  qh_vertexneighbors(qh);
  qh_findgood_all(qh, qh->facet_list);
  qh->RANDOMdist = False;

  facetT *facetlist = qh->facet_list;
  setT *facets = nullptr;
  boolT printall = !qh_ALL;

  // qh_printfacets(fp, qh_PRINToff, facetlist, facets, printall);
  int numcenters, numvertices = 0, numneighbors, numinf, vertex_i, vertex_n;
  facetT *facet, *neighbor, **neighborp;
  setT *vertices;
  vertexT *vertex;
  boolT islower;
  unsigned int numfacets = (unsigned int)qh->num_facets;

  vertices =
      qh_markvoronoi(qh, facetlist, facets, printall, &islower, &numcenters);
  FOREACHvertex_i_(qh, vertices)
  {
    if (vertex) {
      numvertices++;
      numneighbors = numinf = 0;
      FOREACHneighbor_(vertex)
      {
        if (neighbor->visitid == 0)
          numinf = 1;
        else if (neighbor->visitid < numfacets)
          numneighbors++;
      }
      if (numinf && !numneighbors) {
        SETelem_(vertices, vertex_i) = nullptr;
        numvertices--;
      }
    }
  }

  Geometry vcells;
  FORALLfacet_(facetlist)
  {
    if (facet->visitid && facet->visitid < numfacets) {
      if (!facet->normal || !facet->upperdelaunay || !qh->ATinfinity) {
        if (!facet->center)
          facet->center = qh_facetcenter(qh, facet->vertices);
        vcells.add_vert(
            Vec3d(facet->center[0], facet->center[1], facet->center[2]));
      }
      else
        vcells.add_vert(Vec3d(1000, 1000, 1000));
    }
  }
  FOREACHvertex_i_(qh, vertices)
  {
    vector<int> face;
    if (vertex) {
      // qh_order_vertexneighbors(vertex);
      qsort(SETaddr_(vertex->neighbors, vertexT),
            qh_setsize(qh, vertex->neighbors), sizeof(facetT *),
            qh_compare_facetvisit);
    }
    if (vertex) {
      FOREACHneighbor_(vertex)
      {
        if (neighbor->visitid < numfacets)
          face.push_back(neighbor->visitid - 1);
      }
    }
    if (face[0] >= 0)
      vcells.add_face(face);
  }
  qh_settempfree(qh, &vertices);

  qhull_cleanup(qh);
  delete[] points;

  for (unsigned int i = 0; i < vcells.faces().size(); i++) {
    Geometry cell;
    for (unsigned int j = 0; j < vcells.faces(i).size(); j++)
      cell.add_vert(vcells.face_v(i, j));
    cell.add_hull();
    cells->push_back(cell);
  }
  return Status::ok();
}

} // namespace anti
