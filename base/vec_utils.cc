/*
   Copyright (c) 2003, Adrian Rossiter

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
   Name: vec_utils.cc
   Description: vector utilities
   Project: Antiprism - http://www.antiprism.com
*/

#include "math_utils.h"
#include "vec_utils.h"
#include "bbox.h"

vec3d line_plane_intersect(vec3d Q, vec3d n, vec3d P0, vec3d P1, int *where, double eps)
{
   double seg_proj = vdot(n, P1 - P0);
   if(fabs(seg_proj) < eps)
      return vec3d();
   double t = vdot(n, Q - P0) / seg_proj;
   if(where) {
      *where = 0;           // X < P0
      if(t>0) *where=2;     // P0 < X < P1
      if(t>1) *where=4;     // P1 < X
      if(t>0-eps && t<0+eps) *where=1;    // X = P0
      if(t>1-eps && t<1+eps) *where=3;    // X = P1
   }
   return P0 + t*(P1 - P0);
}

vec3d line_plane_intersect(vec3d Q0, vec3d Q1, vec3d Q2, vec3d P0, vec3d P1, int *where, double eps)
{
   vec3d n = vcross(Q1 - Q0, Q2 - Q0);
   //n.unit();
   return line_plane_intersect(Q0, n, P0, P1, where, eps);
}

// http://astronomy.swin.edu.au/~pbourke/geometry/planeplane/
bool two_plane_intersect(vec3d Q0, vec3d n0, vec3d Q1, vec3d n1,
      vec3d &P, vec3d &dir, double eps)
{
   dir = vcross(n0, n1)/ (n0.mag()*n1.mag());
   if(fabs(dir.mag()) < eps)
      return false;
   dir.to_unit();
   double d1 = vdot(Q0, n0);
   double d2 = vdot(Q1, n1);
   double det = vdot(n0, n0)*vdot(n1, n1) - pow(vdot(n0, n1), 2);
   double c1 = (d1*vdot(n1, n1) - d2*vdot(n0, n1)) / det;
   double c2 = (d2*vdot(n0, n0) - d1*vdot(n0, n1)) / det;
   P = c1*n0 + c2*n1;
   return true;
}

// http://astronomy.swin.edu.au/~pbourke/geometry/planes/
bool three_plane_intersect(vec3d Q0, vec3d n0, vec3d Q1, vec3d n1,
      vec3d Q2, vec3d &n2, vec3d &P, double eps)
{
   double tri_prod = vtriple(n0, n1, n2);
   if(fabs(tri_prod)/(n0.mag()*n1.mag()*n2.mag())<eps)
      return false;
   double d1 = vdot(Q0, n0);
   double d2 = vdot(Q1, n1);
   double d3 = vdot(Q2, n2);
   P = (d1*vcross(n1, n2) + d2*vcross(n2, n0) + d3*vcross(n0, n1))/ tri_prod;
   return true;
}

// http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
bool lines_nearest_points(vec3d P0, vec3d P1, vec3d Q0, vec3d Q1,
                          vec3d &P, vec3d &Q, double eps)
{
    vec3d u = (P1-P0);
    vec3d v = (Q1-Q0);
    vec3d w = P0-Q0;
    double a = vdot(u, u);
    double b = vdot(u, v);
    double c = vdot(v, v);
    double d = vdot(u, w);
    double e = vdot(v, w);
    double D = a*c-b*b;
    double s = (b*e-c*d)/D;
    double t = (a*e-b*d)/D;
    P = P0 + s*u;
    Q = Q0 + t*v;
    return (D >= eps);
}


vec3d nearest_point(vec3d P, vec3d Q0, vec3d Q1, vec3d Q2, double eps)
{
   vec3d norm = vcross(Q0 - Q1, Q2 - Q1);
   if(norm.mag2()<=(eps*eps))
      return vec3d();
   
   norm.to_unit();
   return P + vdot(Q0-P, norm)*norm;
}


vec3d nearest_point(vec3d P, const vector <vec3d> &plane, double eps)
{
   switch(plane.size()) {
      case(0):
         return P;
      case(1):
         return plane[0];
      case(2):
         if((plane[1] - plane[0]).mag2() <= eps*eps)
            return plane[0];
         else
            return nearest_point(P, plane[0], plane[1]);
      default:
         vec3d nearest;
         unsigned int p2;
         for(p2=1; p2<plane.size(); p2++)
            if((plane[p2] - plane[0]).mag2() > eps*eps)
               break;
         if(p2==plane.size())
             return plane[0];
         if(p2==plane.size()-1)
            return nearest_point(P, plane[0], plane[1]);
         
         for(unsigned int i=p2; i<plane.size(); i++) {
            nearest = nearest_point(P, plane[0], plane[p2], plane[i], eps);
            if(nearest.is_set())
               return nearest;
         }
         
         return nearest_point(P, plane[0], plane[1]);
   }         
}

vec3d nearest_point(vec3d P, const vector<vec3d> &plane, const vector<int> &idxs, double eps)
{
   vector<vec3d> plane2;
   for(unsigned int i=0; i< idxs.size(); i++)
      plane2.push_back(plane[idxs[i]]);

   return nearest_point(P, plane2, eps);
}


double angle_around_axis(const vec3d &v0, const vec3d &v1, const vec3d &axis)
{
   vec3d n0 = vcross(v0, axis).unit();
   vec3d n1 = vcross(v1, axis).unit();
   double ang = acos(safe_for_trig(vdot(n0, n1)));
   if(vdot(axis, vcross(n0, n1))<0)
      ang = 2*M_PI - ang;
   return ang;
}


vec3d lines_intersection(const vec3d &P0, const vec3d &P1, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d P, Q;
   if(eps==0.0) {                  // always return an intersection point
      if(!lines_nearest_points(P0, P1, Q0, Q1, P, Q)) {
         P = (P0+P1)/2.0;
         Q = (Q0+Q1)/2.0;
      }
      return (P+Q)/2.0;
   }

   // lines might not be parallel and still miss so check if nearest points is not zero
   if(lines_nearest_points(P0, P1, Q0, Q1, P, Q, eps) && ((P-Q).mag() < eps))
      return (P+Q)/2.0;
   else
      return vec3d();
}


bool in_segment(const vec3d &P, const vec3d &Q0, const vec3d &Q1, double eps)
{
   bound_box bb;
   vector<vec3d> points(2);
   
   points[0] = Q0;
   points[1] = Q1;

   bb.add_points(points);
   vec3d min = bb.get_min();
   vec3d max = bb.get_max();

   return (compare(min,P,eps)<=0 && compare(max,P,eps)>=0);
}

vec3d lines_intersection_in_segments(const vec3d &P0, const vec3d &P1, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d intersection = lines_intersection(P0, P1, Q0, Q1, eps);
   
   if (intersection.is_set() && (!in_segment(intersection, P0, P1, eps) || !in_segment(intersection, Q0, Q1, eps)))
      intersection = vec3d();

   return intersection;
}

vec3d point_in_segment(const vec3d &P, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d intersection = nearest_point(P, Q0, Q1);

   if (intersection.is_set() && (compare(P,intersection,eps) || !in_segment(intersection, Q0, Q1, eps)))
      intersection = vec3d();

   return intersection;
   //return (intersection.is_set() && !compare(P,intersection,eps) && in_segment(intersection, Q0, Q1, eps));
}

