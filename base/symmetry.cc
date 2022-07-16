/*
   Copyright (c) 2008-2021, Adrian Rossiter

   Project: Antiprism - http://www.antiprism.com

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

/* \file symmetry.cc
 *\brief symmetry handling
 */

#include "symmetry.h"
#include "geometryinfo.h"
#include "mathutils.h"
#include "utils.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>
#include <numeric>
#include <set>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::swap;
using std::vector;

namespace anti {

const Vec3d A3(1, 1, 1); // 3-fold axis for T, O and I symmetry types
const Vec3d A5(0, 1, (sqrt(5) + 1) / 2); // 5-fold axis for I symmetry type

Transformations &Transformations::product(const Transformations &s1,
                                          const Transformations &s2)
{
  clear();
  set<Trans3d>::iterator i1, i2;
  for (i1 = s1.trans.begin(); i1 != s1.trans.end(); i1++)
    for (i2 = s2.trans.begin(); i2 != s2.trans.end(); i2++)
      add((*i1) * (*i2));
  return *this;
}

Transformations &Transformations::product_with(const Transformations &s)
{
  SymTransSet trans_orig = trans;
  clear();
  set<Trans3d>::iterator i0, i1;
  for (i0 = trans_orig.begin(); i0 != trans_orig.end(); i0++)
    for (i1 = s.trans.begin(); i1 != s.trans.end(); i1++)
      add((*i0) * (*i1));
  return *this;
}

Transformations &Transformations::conjugate(const Trans3d &t)
{
  SymTransSet conj;
  Trans3d inv = t.inverse();
  for (const auto &tran : trans)
    conj.insert(t * tran * inv);
  trans = conj;
  return *this;
}

struct SymTransCmp {
  bool operator()(const Trans3d &t1, const Trans3d &t2)
  {
    return compare(t1, t2, sym_eps);
  }
};
Transformations &Transformations::intersection(const Transformations &s1,
                                               const Transformations &s2)
{
  clear();
  set_intersection(s1.trans.begin(), s1.trans.end(), s2.trans.begin(),
                   s2.trans.end(), inserter(trans, trans.end()),
                   SymTransLess());
  return *this;
}

Transformations &Transformations::subtract(const Transformations &s)
{
  SymTransSet diff;
  set_difference(trans.begin(), trans.end(), s.trans.begin(), s.trans.end(),
                 inserter(diff, diff.begin()), SymTransLess());
  trans = diff;
  return *this;
}

Transformations Transformations::lcoset(const Trans3d &t) const
{
  Transformations coset;
  for (const auto &tran : trans)
    coset.add(t * tran);
  return coset;
}

Transformations operator*(const Transformations &s1, const Transformations &s2)
{
  return Transformations().product(s1, s2);
}

Transformations &operator*=(Transformations &s1, const Transformations &s2)
{
  return s1.product_with(s2);
}

int Transformations::lcosets(Transformations sub,
                             vector<Transformations> &lcosets) const
{
  lcosets.clear();
  Transformations whole = *this;
  sub += Trans3d(); // must include unit!
  while (whole.is_set()) {
    Trans3d tr = *whole.trans.begin(); // select a transf from init list
    lcosets.push_back(sub.lcoset(tr)); // find equivalent transformations
    whole.subtract(lcosets.back());    // remove them
  }
  return lcosets.size();
}

Transformations operator+(const Transformations &s, const Trans3d &m)
{
  return Transformations(s).add(m);
}

Transformations &operator+=(Transformations &s, const Trans3d &m)
{
  return s.add(m);
}

int compare(const Transformations &t0, const Transformations &t1)
{
  if (t0.size() > t1.size())
    return 1;
  else if (t0.size() < t1.size())
    return -1;

  // same size so test matrices one at a time
  auto i0 = t0.begin();
  auto i1 = t1.begin();
  for (; i0 != t0.end(); ++i0, ++i1) {
    int cmp = compare(*i0, *i1, sym_eps);
    if (cmp)
      return cmp;
  }

  return 0;
}

Transformations &Transformations::min_set(const Transformations &tr_whole,
                                          const Transformations &tr_part,
                                          const Trans3d &pos)
{
  Transformations whole = tr_whole;
  Transformations part = tr_part;
  part.conjugate(pos);
  Transformations inter;
  inter.intersection(part, whole);
  inter += Trans3d(); // must include unit!

  while (whole.is_set()) {
    Trans3d tr = *whole.trans.begin();        // select a transf from init list
    add(tr);                                  // and add it to the final list
    Transformations coset = inter.lcoset(tr); // find equivalent transformations
    whole.subtract(coset);                    // remove them
  }
  return *this;
}

/// Class for working with Schoenflies notation
class sch_gen : public Transformations {
private:
  sch_gen &operator*(const Transformations &s)
  {
    product_with(s);
    return *this;
  }
  sch_gen &operator+(const Trans3d &m)
  {
    add(m);
    return *this;
  }

  /// Set up mirror transformation group.
  /**\param norm The normal for the mirror.
   *\return reference to this object with the transformations set. */
  sch_gen &refl(const Vec3d &norm)
  {
    return unit() + Trans3d::reflection(norm);
  }

public:
  /// Constructor
  /**\param t initialise with these transformations. */
  sch_gen(Transformations t = Transformations()) : Transformations(t){};

  /// Set up unit symmetry transformation group.
  /**No relevant alignment.
   *\return reference to this object with the transformation set. */
  sch_gen &unit()
  {
    add(Trans3d());
    return *this;
  }

  /// Set up horizontal mirror transformation group.
  /**Mirror normal (0,0,1).
   *\return reference to this object with the transformations set. */
  sch_gen &h_refl() { return refl(Vec3d::Z); }

  /// Set up vertical mirror transformation group.
  /**Mirror normal (0,1,0).
   *\return reference to this object with the transformations set. */
  sch_gen &v_refl() { return refl(Vec3d::Y); }

  /// Set up vertical mirror transformation group.
  /**Mirror normal (0,1,0) rotated PI/2n radians around (0,0,1).
   *\param n used for angle to rotate mirror, as PI/n radians.
   *\return reference to this object with the transformations set. */
  sch_gen &v_refl(int n)
  {
    return refl(Trans3d::rotate(Vec3d::Z, -0.5 * M_PI / n) * Vec3d::Y);
  }

  /// Set up dihedral symmetry transformation group.
  /**Rotation axis (1,0,0).
   *\return reference to this object with the transformations set. */
  sch_gen &C2()
  {
    C(2);
    conjugate(Trans3d::rotate(Vec3d::Z, Vec3d::X));
    return *this;
  }

  /// Set up Cs transformation group.
  /**Mirror normal in direction (0,0,1).
   *\return reference to this object with the transformations set. */
  sch_gen &Cs() { return h_refl(); }

  /// Set up Ci transformation group.
  /**No relevant alignment.
   *\return reference to this object with the transformations set. */
  sch_gen &Ci() { return unit() + Trans3d::inversion(); }

  /// Set up C symmetry transformation group.
  /**Principal axis (0,0,1).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &C(int n);

  /// Set up Cv symmetry transformation group.
  /**Principal axis (0,0,1). Mirror normal (0,1,0).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &Cv(int n) { return C(n) * sch_gen().v_refl(); }

  /// Set up Ch symmetry transformation group.
  /**Principal axis (0,0,1). Mirror normal in direction (0,0,1).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &Ch(int n) { return C(n) * sch_gen().h_refl(); }

  /// Set up D symmetry transformation group.
  /**Principal axis (0,0,1). Dihedral axis in direction (1,0,0).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &D(int n) { return C(n) * sch_gen().C2(); }

  /// Set up Dv symmetry transformation group.
  /**Principal axis (0,0,1). Dihedral axis (1,0,0).
   * Vertical mirror normal (0,1,0) rotated PI/2n radians around (0,0,1).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &Dv(int n) { return D(n) * sch_gen().v_refl(n); }

  /// Set up Dh symmetry transformation group.
  /**Principal axis (0,0,1). Dihedral axis (1,0,0).
   * Horizontal mirror normal (0,0,1).
   *\param n principal axis is n-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &Dh(int n) { return D(n) * sch_gen().h_refl(); }

  /// Set up Sn symmetry transformation group.
  /**Principal axis (0,0,1). Horizontal mirror normal (0,0,1).
   *\param n principal axis is rotational n/2-fold.
   *\return reference to this object with the transformations set. */
  sch_gen &S(int n)
  {
    return C(n / 2) *
           (sch_gen().unit() + Trans3d::reflection(Vec3d::Z) *
                                   Trans3d::rotate(Vec3d::Z, 2 * M_PI / n));
  }

  /// Set up T symmetry transformation group.
  /**3-fold axes (1,1,1), (1,-1,-1).
   *\return reference to this object with the transformations set. */
  sch_gen &T()
  {
    return D(2) * (sch_gen().C(3).conjugate(Trans3d::rotate(Vec3d::Z, A3)));
  }

  /// Set up Td symmetry transformation group.
  /**3-fold axes (1,1,1), (1,-1,-1).
   *\return reference to this object with the transformations set. */
  sch_gen &Td() { return T() * sch_gen().refl(Vec3d::Y + Vec3d::X); }

  /// Set up Th symmetry transformation group.
  /**3-fold axes (1,1,1), (1,-1,-1).
   *\return reference to this object with the transformations set. */
  sch_gen &Th() { return T() * sch_gen().Ci(); }

  /// Set up O symmetry transformation group.
  /**4-fold axes (1,0,0), (0,1,0).
   *\return reference to this object with the transformations set. */
  sch_gen &O()
  {
    return T() * (sch_gen().unit() + Trans3d::rotate(Vec3d::X, M_PI / 2));
  }

  /// Set up Oh symmetry transformation group.
  /**4-fold axes (1,0,0), (0,1,0).
   *\return reference to this object with the transformations set. */
  sch_gen &Oh()
  {
    return O() * sch_gen().Ci();
    return *this;
  }

  /// Set up I symmetry transformation group.
  /**5-fold axes (0,1,phi), (0,1,-phi).
   *\return reference to this object with the transformations set. */
  sch_gen &I()
  {
    return T() * sch_gen().C(5).conjugate(Trans3d::rotate(Vec3d::Z, A5));
  }

  /// Set up Ih symmetry transformation group.
  /**5-fold axes (0,1,phi), (0,1,-phi).
   *\return reference to this object with the transformations set. */
  sch_gen &Ih() { return I() * sch_gen().Ci(); }
};

sch_gen &sch_gen::C(int n)
{
  clear();
  unit();
  for (int i = 1; i < n; i++)
    add(Trans3d::rotate(Vec3d::Z, 2 * M_PI * i / n));
  return *this;
}

void Isometry::normalise()
{
  if (axis.is_set()) {
    if ((axis[0] < -sym_eps) || (axis[0] < sym_eps && axis[1] < -sym_eps) ||
        (axis[0] < sym_eps && axis[1] < sym_eps && axis[2] < -sym_eps)) {
      axis = -axis;
      ang = -ang;
    }
    if (fabs(axis[0]) > sym_eps || fabs(axis[1]) > sym_eps ||
        fabs(axis[2]) > sym_eps)
      axis.to_unit();
    else
      axis = Vec3d();
  }

  ang = fmod(ang, 2 * M_PI);
  if (ang < 0) {
    if (ang > -sym_eps)
      ang = 0;
    else
      ang += 2 * M_PI;
  }
}

bool Isometry::is_isometry(Trans3d m) const
{
  Vec3d rows[3];
  rows[0] = Vec3d(m[0], m[1], m[2]);
  rows[1] = Vec3d(m[4], m[5], m[6]);
  rows[2] = Vec3d(m[8], m[9], m[10]);
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++) {
      double dot = vdot(rows[i], rows[j]);
      /// dot should equal 1 if i==j, otherwise 0
      if (dot < (i == j) - sym_eps || dot > (i == j) + sym_eps)
        return false;
    }

  return true;
}

Isometry &Isometry::init(Trans3d m)
{
  if (!is_isometry(m)) {
    rot_type = rt_none;
    return *this;
  }

  transl = Vec3d(m[3], m[7], m[11]); // translation column
  m[3] = m[7] = m[11] = 0;           // zero translation column
  double det = m.det();
  if (det < 0)
    m *= Trans3d::inversion();

  Vec4d quat = m.get_quaternion();
  double cos_a = quat[3];
  ang = 2 * acos(safe_for_trig(cos_a));

  if (sqrt(1 - cos_a * cos_a) > anti::epsilon)
    axis = Vec3d(quat[0], quat[1], quat[2]);
  else
    axis.unset();

  if (det > 0) {
    if (fabs(ang) < sym_eps)
      rot_type = rt_unit;
    else
      rot_type = rt_rot;
  }
  else {
    if (fabs(ang) < sym_eps)
      rot_type = rt_inv;
    else if (fabs(fabs(ang) - M_PI) < sym_eps)
      rot_type = rt_refl;
    else
      rot_type = rt_rot_refl;
    ang = M_PI + ang;
  }

  normalise();
  return *this;
}

void Isometry::dump() const
{
  const char *et[] = {"none", "unit", "rot", "inv", "refl", "rot_refl"};
  fprintf(stderr, "rot_type=%s", et[rot_type]);
  if (rot_type == rt_rot || rot_type == rt_refl || rot_type == rt_rot_refl) {
    fprintf(stderr, ", axis=(%-.3f, %-.3f, %-.3f)", axis[0], axis[1], axis[2]);
    if (rot_type == rt_rot || rot_type == rt_rot_refl)
      fprintf(stderr, ", ang=%3.0f", ang * 180 / M_PI);
  }
  if (rot_type != rt_none)
    fprintf(stderr, ", transl=(%-.3f, %-.3f, %-.3f)", transl[0], transl[1],
            transl[2]);
  fprintf(stderr, "\n");
}

SymmetryAxis::SymmetryAxis(const Trans3d &m)
{
  Isometry rot(m);
  double ang = rot.get_ang();
  if (fabs(ang) < sym_eps)
    nfold = 1;
  else {
    long tmp;
    double2rational(ang / (2 * M_PI), tmp, nfold, sym_eps);
  }
  axis = rot.get_axis();
  // map index {rt_none=0, rt_unit, rt_rot, rt_inv, rt_refl, rt_rot_refl}
  int sym_map[] = {Symmetry::C1, Symmetry::C1, Symmetry::C,
                   Symmetry::Ci, Symmetry::Cs, Symmetry::S};
  sym_type = sym_map[rot.get_rot_type()];
}

void SymmetryAxis::dump() const
{
  fprintf(stderr, "axis=%s", Symmetry(sym_type, nfold).get_symbol().c_str());
  if (axis.is_set())
    fprintf(stderr, ": a=(% .3f,% .3f,% .3f)", axis[0], axis[1], axis[2]);
  if (perp.is_set())
    fprintf(stderr, ", p=(% .3f,% .3f,% .3f)", perp[0], perp[1], perp[2]);
  fprintf(stderr, "\n");
}

bool SymmetryAxis::operator<(const SymmetryAxis &s) const
{
  int cmp = compare(axis, s.axis, sym_eps);
  if (cmp < 0)
    return true;
  else if (cmp == 0) {
    if (sym_type < s.sym_type)
      return true;
    else if (sym_type == s.sym_type) {
      if (nfold < s.nfold)
        return true;
      else if (nfold == s.nfold) {
        int cmp = compare(perp, s.perp, sym_eps);
        if (cmp < 0)
          return true;
      }
    }
  }
  return false;
}

//-----------------------------------------------------------------
// Subspace

Subspace::Subspace(SubspaceType type, Vec3d point, Vec3d dir)
    : type(type), point(point)
{
  if (dir.is_set())
    direction = dir.unit();
}

Vec3d Subspace::nearest_point(Vec3d P) const
{
  if (type == SubspaceType::none)
    return Vec3d();
  else if (type == SubspaceType::point)
    return point;
  else if (type == SubspaceType::line)
    return anti::nearest_point(P, point, point + direction);
  else if (type == SubspaceType::plane)
    return P + vdot(point - P, direction) * direction;
  else // type == SubspaceType::space
    return P;
}

//-----------------------------------------------------------------
// Symmetry

// find symmetry from rotational axes
void Symmetry::find_full_sym_type(const set<SymmetryAxis> &full_sym)
{
  // find two greatest n-fold axes, do not allow axes to be colinear
  SymmetryAxis max_fold1, max_fold2;
  bool has_dv = false, has_dh = false;
  for (const auto &si : full_sym) {
    int nfold = si.get_nfold();
    if (nfold > max_fold1.get_nfold())
      max_fold1 = si;
    else if (nfold > max_fold2.get_nfold() &&
             compare(si.get_axis(), -max_fold1.get_axis(), sym_eps) != 0)
      max_fold2 = si;
    if (si.get_sym_type() == Symmetry::Dv)
      has_dv = true;
    if (si.get_sym_type() == Symmetry::Dh)
      has_dh = true;
  }
  nfold = max_fold1.get_nfold();

  // dihedral axes only
  if (max_fold1.get_nfold() == 2) {
    SymmetryAxis ax = max_fold1;
    for (const auto &si : full_sym)
      if (si.get_sym_type() == Symmetry::Dv) {
        ax = si;
        break;
      }
    sym_type = ax.get_sym_type();
    to_std = Trans3d::align(ax.get_axis(), ax.get_perp(), Vec3d::Z, Vec3d::X);
  }

  // principal axis
  else if (max_fold1.get_nfold() > 2 && max_fold2.get_nfold() <= 2) {
    sym_type = max_fold1.get_sym_type();
    Vec3d perp = max_fold1.get_perp();
    if (perp.is_set())
      to_std = Trans3d::align(max_fold1.get_axis(), max_fold1.get_perp(),
                              Vec3d::Z, Vec3d::X);
    else
      to_std = Trans3d::rotate(max_fold1.get_axis(), Vec3d::Z);
  }

  // tetrahedral
  else if ((max_fold1.get_nfold() == 3 || max_fold1.get_nfold() == 6) &&
           // axes may be part of incompletely detected O or I symmetry,
           // in which case leave to be detected as C1
           double_eq(fabs(vdot(max_fold1.get_axis(), max_fold2.get_axis())),
                     1.0 / 3, sym_eps)) {
    if (has_dh)
      sym_type = Symmetry::Th;
    else if (has_dv)
      sym_type = Symmetry::Td;
    else
      sym_type = Symmetry::T;

    Vec3d A3b = Trans3d::rotate(Vec3d::Z, M_PI) * A3;
    if (vdot(max_fold1.get_axis(), max_fold2.get_axis()) > 0)
      A3b *= -1;

    to_std =
        Trans3d::align(max_fold1.get_axis(), max_fold2.get_axis(), A3, A3b);
  }

  // octahedral
  else if (max_fold1.get_nfold() == 4) {
    if (has_dh)
      sym_type = Symmetry::Oh;
    else
      sym_type = Symmetry::O;

    to_std = Trans3d::align(max_fold1.get_axis(), max_fold2.get_axis(),
                            Vec3d::Z, Vec3d::X);
  }

  // icosahedral
  else if (max_fold1.get_nfold() == 5) {
    if (has_dh)
      sym_type = Symmetry::Ih;
    else
      sym_type = Symmetry::I;

    Vec3d A5b = Trans3d::rotate(A3, 2 * M_PI / 3) * A5;
    if (vdot(max_fold1.get_axis(), max_fold2.get_axis()) < 0)
      A5b *= -1;

    to_std =
        Trans3d::align(max_fold1.get_axis(), max_fold2.get_axis(), A5, A5b);
  }

  // unknown symmetry
  else {
    // sym_type = Symmetry::unknown;
    // Don't allow to fail
    axes.clear(); // clear any detected axes
    sym_type = Symmetry::C1;
    to_std = Trans3d();
  }
}

static bool operator<(const Vec3d &v1, const Vec3d &v2)
{
  return compare(v1, v2, sym_eps) < 0;
}

Symmetry::Symmetry(const Transformations &ts)
    : sym_type(unknown), nfold(1), to_std(Trans3d())
{
  if (!ts.size())
    return;

  // Translate the set of symmetry elements to fix the origin
  Vec3d fixed_pt(0, 0, 0);
  for (const auto &t : ts)
    fixed_pt += Vec3d(t[3], t[7], t[11]);
  fixed_pt /= ts.size(); // centroid of points where origin is sent

  Trans3d transl = Trans3d::translate(-fixed_pt);
  Transformations o_ts = ts;
  o_ts.conjugate(transl);

  for (const auto &o_t : o_ts)
    axes.insert(SymmetryAxis(o_t));

  bool inv = false;
  vector<Vec3d> refls, dihs;
  for (const auto &axe : axes) {
    if (axe.get_sym_type() == Symmetry::Cs)
      refls.push_back(axe.get_axis());
    else if (axe.get_sym_type() == Symmetry::C && axe.get_nfold() == 2)
      dihs.push_back(axe.get_axis());
    else if (axe.get_sym_type() == Symmetry::Ci)
      inv = true;
  }

  map<Vec3d, vector<SymmetryAxis>> v2ax;
  for (const auto &axe : axes) {
    const Vec3d &axis = axe.get_axis();
    if (axis.is_set())
      v2ax[axis].push_back(axe);
  }

  if (v2ax.size() == 0) { // no axes
    axes.clear();
    if (inv)
      sym_type = Symmetry::Ci;
    else
      sym_type = Symmetry::C1;
  }
  else if (v2ax.size() == 1 && v2ax.begin()->second.size() == 1 &&
           v2ax.begin()->second[0].get_sym_type() == Symmetry::Cs) { // one refl
    axes.clear();
    to_std = Trans3d::rotate(v2ax.begin()->second[0].get_axis(), Vec3d::Z);
    sym_type = Symmetry::Cs;
  }
  else { // remaining possibilities have rotational axis
    set<SymmetryAxis> full_sym;
    map<Vec3d, vector<SymmetryAxis>>::iterator mi;
    for (mi = v2ax.begin(); mi != v2ax.end(); mi++) {
      Vec3d axis = mi->first;
      int nfold = 0;
      int sym = -1;
      int h_refl = false;
      for (auto &vi : mi->second) {
        if (vi.get_sym_type() == Symmetry::Cs)
          h_refl = true;
        if (vi.get_nfold() > nfold) {
          nfold = vi.get_nfold();
          sym = vi.get_sym_type();
        }
      }

      SymmetryAxis sym_ax;
      sym_ax.set_axis(axis);
      if (!h_refl) { // no horz reflection -> S or C
        sym_ax.set_sym_type(sym);
        sym_ax.set_nfold(nfold);
      }
      else { // horz refl -> Ch
        if (sym_ax.get_sym_type() == Symmetry::S)
          sym_ax.set_nfold(nfold / 2);
        else
          sym_ax.set_nfold(nfold);
        sym_ax.set_sym_type(Symmetry::Ch);
      }

      // vertical reflection

      Trans3d align_axis = Trans3d::rotate(axis, Vec3d::Z);
      if (sym_ax.get_sym_type() == Symmetry::S)
        align_axis =
            Trans3d::rotate(Vec3d::Z, M_PI / sym_ax.get_nfold()) * align_axis;
      Vec3d closest_to_y;
      for (auto &refl : refls) {
        if (fabs(vdot(refl, axis)) < sym_eps) {
          Vec3d v = align_axis * refl;
          if (!closest_to_y.is_set() || fabs(v[1]) > closest_to_y[1])
            closest_to_y = (1 - 2 * (v[1] < 0)) * refl;
        }
      }
      if (closest_to_y.is_set()) {
        if (sym_ax.get_sym_type() == Symmetry::C) {
          sym_ax.set_sym_type(Symmetry::Cv);
          sym_ax.set_perp(Trans3d::rotate(axis, -M_PI / 2) * closest_to_y);
        }
        else if (sym_ax.get_sym_type() == Symmetry::Ch) {
          sym_ax.set_sym_type(Symmetry::Dh);
          // sym_ax.set_perp(Trans3d::rotate(axis, M_PI/2+M_PI/nfold)*(*vi));
          sym_ax.set_perp(Trans3d::rotate(axis, -M_PI / 2) * closest_to_y);
        }
        else { // Symmetry::S
          sym_ax.set_sym_type(Symmetry::Dv);
          sym_ax.set_nfold(nfold / 2);
          sym_ax.set_perp(
              Trans3d::rotate(axis,
                              M_PI / 2 - 0.5 * M_PI / sym_ax.get_nfold()) *
              closest_to_y);
        }
      }

      // dihedral axis
      if (sym_ax.get_sym_type() != Symmetry::Dv &&
          sym_ax.get_sym_type() != Symmetry::Dh) {
        for (auto &dih : dihs) {
          if (fabs(vdot(dih, axis)) < sym_eps) {
            if (sym_ax.get_sym_type() == Symmetry::C) {
              sym_ax.set_perp(dih);
              sym_ax.set_sym_type(Symmetry::D);
            }
            else { // sym==Symmetry::Cv
              sym_ax.set_sym_type(Symmetry::Dv);
              sym_ax.set_perp(
                  Trans3d::rotate(axis, -0.5 * M_PI / sym_ax.get_nfold()) *
                  (sym_ax.get_perp()));
            }
            break;
          }
        }
      }
      if (nfold > 1) { // has rotation axis (not just refl)
        full_sym.insert(sym_ax);
        // sym_ax.dump();
      }
    }

    axes = full_sym;
    find_full_sym_type(full_sym);
  }

  to_std *= transl; // first move the fixed point to the origin;
}

static inline int edge_seen(map<pair<int, int>, unsigned char> &e_seen, int v0,
                            int v1, bool update)
{
  unsigned char dir = 1;
  if (v0 > v1) {
    swap(v0, v1);
    dir = 2;
  }
  pair<int, int> edge(v0, v1);
  auto ei = e_seen.find(edge);
  if (ei == e_seen.end()) { // not traversed in either direction
    if (update)
      e_seen[edge] = dir; // mark as traversed in this direction
    return 0;             // 0: unseen
  }
  if (ei->second & dir)
    return 2; // 2: seen and already traversed in this direction
  if (update)
    ei->second |= dir; // mark as traversed in this direction
  return 1;            // 1: seen but not traversed in this direction
}

static inline int edge_check(map<pair<int, int>, unsigned char> &e_seen, int v0,
                             int v1)
{
  return edge_seen(e_seen, v0, v1, false);
}

static inline int edge_mark(map<pair<int, int>, unsigned char> &e_seen, int v0,
                            int v1)
{
  return edge_seen(e_seen, v0, v1, true);
}

// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.30.6536
// Symmetries of Polyhedra: Detection and Applications
// by X. Y. Jiang, H. Bunke
// ftp://ftp.iam.unibe.ch/pub/TechReports/1994/iam-94-012.ps.gz
// 3.1.1 The algorithm of Jiang & Bunke

static int find_path(vector<int> &path, vector<int> &v_code,
                     const vector<int> &edge, const vector<vector<int>> &v_cons,
                     const vector<int> *test_path = nullptr,
                     const vector<int> *test_v_code = nullptr)
{
  path.clear();
  v_code = vector<int>(v_cons.size(), -1);
  map<pair<int, int>, unsigned char> e_seen;
  int v_cnt = 0;
  int v_cur = edge[0];
  path.push_back(v_cur);
  v_code[v_cur] = v_cnt++;
  int v_next = edge[1];
  while (true) {
    int seen = edge_mark(e_seen, v_cur, v_next);
    int v_new = -1;
    vector<int>::const_iterator vi, vi_orig;
    if (v_code[v_next] < 0) { // new vertex, exit from "right"
      vi = find(v_cons[v_next].begin(), v_cons[v_next].end(), v_cur);
      if (++vi == v_cons[v_next].end())
        vi = v_cons[v_next].begin();
      v_code[v_next] = v_cnt++;
      v_new = *vi;
    }
    else {        // previously seen vertex
      if (seen) { // previously seen edge, exit from "right"
        vi_orig = find(v_cons[v_next].begin(), v_cons[v_next].end(), v_cur);
        for (vi = vi_orig + 1; vi != vi_orig; ++vi) { // check for an edge out
          if (vi == v_cons[v_next].end()) {
            vi = v_cons[v_next].begin(); // wrap to beginning
            if (vi == vi_orig)
              break;
          }
          if (edge_check(e_seen, v_next, *vi) < 2) {
            v_new = *vi;
            break;
          }
        }
        if (vi == vi_orig) // finished, didn't find an edge to leave from
          return 1;
      }
      else { // new edge, go back
        v_new = v_cur;
      }
    }

    // stop if the coded path is different to the test coded path
    if (test_path &&
        v_code[v_next] != (*test_v_code)[(*test_path)[path.size()]])
      return 0;

    path.push_back(v_next);
    v_cur = v_next;
    v_next = v_new;
  }
}

static void update_equiv_elems(vector<map<int, set<int>>> &equiv_elems,
                               const vector<map<int, set<int>>> &new_equivs,
                               int cnts[3])
{
  for (int i = 0; i < 3; i++) {
    map<int, set<int>>::const_iterator mi_new;
    for (mi_new = new_equivs[i].begin(); mi_new != new_equivs[i].end();
         ++mi_new) {
      int to = *mi_new->second.begin();
      if (to >= cnts[i])
        to -= cnts[i];
      for (int si : mi_new->second) {
        int from = (si < cnts[i]) ? si : si - cnts[i];
        equiv_elems[i][to].insert(from);
      }
    }
  }
}

static void equiv_elems_to_sets(vector<vector<set<int>>> &equiv_sets,
                                vector<map<int, set<int>>> &equiv_elems,
                                vector<map<int, set<int>>> &orig_equivs)
{
  equiv_sets.clear();
  equiv_sets.resize(3);
  for (int i = 0; i < 3; i++) {
    // map each index to the to minimum index it is equivalent to
    map<int, int> idx_to_min_idx;
    for (auto &mi : equiv_elems[i]) {
      auto si = mi.second.begin();
      int min_idx = *si;
      for (; si != mi.second.end(); ++si) {
        auto idx_i = idx_to_min_idx.find(*si);
        if (idx_i == idx_to_min_idx.end())
          idx_to_min_idx[*si] = min_idx;
        else if (idx_i->second > min_idx)
          idx_i->second = min_idx;
      }
    }

    // map each minimum index to its set of equivalent indexes
    map<int, set<int>> min_idx_to_idxs;
    map<int, int>::const_iterator idx_i;
    for (idx_i = idx_to_min_idx.begin(); idx_i != idx_to_min_idx.end(); ++idx_i)
      min_idx_to_idxs[idx_i->second].insert(idx_i->first);

    // Copy sets to the return argument
    map<int, set<int>>::const_iterator midx_i;
    for (midx_i = min_idx_to_idxs.begin(); midx_i != min_idx_to_idxs.end();
         ++midx_i) {
      equiv_sets[i].push_back(set<int>());
      for (int si : midx_i->second)
        for (auto si2 = equiv_elems[i][si].begin();
             si2 != equiv_elems[i][si].end(); ++si2)
          equiv_sets[i].back().insert(orig_equivs[i][*si2].begin(),
                                      orig_equivs[i][*si2].end());
    }
    std::sort(equiv_sets[i].begin(), equiv_sets[i].end());
  }
}

static bool is_sym(const Geometry &test_geom, const Geometry &geom,
                   const vector<int> &test_v_code, const vector<int> &v_code,
                   bool orient, Trans3d &trans,
                   vector<map<int, set<int>>> &new_equivs)
{
  int v_sz = test_geom.verts().size();
  // code to vertex idx for this sym
  vector<int> c2v_map(v_sz);
  for (int v = 0; v < v_sz; v++)
    c2v_map[v_code[v]] = v;
  vector<int> v_map(v_sz);
  for (int v = 0; v < v_sz; v++)
    v_map[v] = c2v_map[test_v_code[v]];

  vector<Vec3d> t_pts(3), pts(3);
  for (int i = 0; i < 2; i++) {
    t_pts[i] = test_geom.verts(i);
    pts[i] = test_geom.verts(v_map[i]);
  }
  // Choose a third point that is not colinear with the first two
  for (int i = 2; i < v_sz; i++) {
    t_pts[2] = test_geom.verts(i);
    pts[2] = test_geom.verts(v_map[i]);
    Vec3d norm = vcross(pts[1] - pts[0], pts[2] - pts[0]);
    if (norm.len2() > (anti::epsilon * anti::epsilon))
      break;
  }

  if (orient)
    transform(pts, Trans3d::inversion());
  trans = Trans3d::align(t_pts, pts);
  if (orient)
    trans = Trans3d::inversion() * trans;
  Geometry s_geom = geom;
  s_geom.transform(trans);

  return check_coincidence(geom, s_geom, &new_equivs, sym_eps);
}

static void set_equiv_elems_identity(const Geometry &geom,
                                     vector<vector<set<int>>> *equiv_sets)
{
  int cnts[3] = {(int)geom.verts().size(), (int)geom.edges().size(),
                 (int)geom.faces().size()};
  equiv_sets->clear();
  equiv_sets->resize(3);
  for (int i = 0; i < 3; i++) {
    (*equiv_sets)[i].resize(cnts[i]);
    for (int j = 0; j < cnts[i]; j++)
      (*equiv_sets)[i][j].insert(j);
  }
}

static int find_syms(const Geometry &geom, Transformations &ts,
                     vector<vector<set<int>>> *equiv_sets)
{
  ts.clear();

  Geometry merged_geom = geom;
  vector<map<int, set<int>>> orig_equivs;
  merge_coincident_elements(merged_geom, "vef", &orig_equivs, anti::epsilon);
  Geometry test_geom = merged_geom;

  GeometryInfo inf(merged_geom);
  // if(inf.num_parts()>1 || !inf.is_orientable()) // for octahemioctaheron=Td
  int dim;
  test_geom.set_hull(msg_str("-A%.15f", 1.0 - sym_eps), &dim);
  if (dim < 2) // contains an infinite axis, can't currently handle this
    return 0;
  test_geom.orient();

  GeometryInfo g_inf(test_geom);
  const vector<vector<int>> &edges = g_inf.get_impl_edges();

  vector<vector<int>> v_cons(g_inf.num_verts());
  for (unsigned int i = 0; i < v_cons.size(); i++) {
    if (g_inf.get_vert_figs()[i].size()) // shouldn't fail for 3d hull
      v_cons[i] = g_inf.get_vert_figs()[i][0];
    else
      v_cons[i] = g_inf.get_vert_cons()[i]; // needed for 2d hull
  }

  vector<vector<int>> r_cons = v_cons;
  for (auto &r_con : r_cons)
    reverse(r_con.begin(), r_con.end());
  const vector<vector<int>> *cons[] = {&v_cons, &r_cons};

  vector<map<int, set<int>>> equiv_elems(3);
  int cnts[3] = {(int)merged_geom.verts().size(),
                 (int)merged_geom.edges().size(),
                 (int)merged_geom.faces().size()};
  vector<int> test_path, path;
  vector<int> test_v_code, v_code;
  find_path(test_path, test_v_code, *edges.begin(), v_cons);
  vector<vector<int>>::const_iterator ei;
  for (ei = edges.begin(); ei != edges.end(); ei++) {
    vector<int> edge = *ei;
    for (int i = 0; i < 2; i++) {
      if (i)
        swap(edge[0], edge[1]);
      for (int orient = 0; orient < 2; orient++) {
        if (find_path(path, v_code, edge, *cons[orient], &test_path,
                      &test_v_code)) {
          Trans3d trans;
          vector<map<int, set<int>>> new_equivs;
          if (is_sym(test_geom, merged_geom, test_v_code, v_code, orient, trans,
                     new_equivs)) {
            ts.add(trans);
            if (equiv_sets)
              update_equiv_elems(equiv_elems, new_equivs, cnts);
          }
        }
      }
    }
  }

  if (equiv_sets)
    equiv_elems_to_sets(*equiv_sets, equiv_elems, orig_equivs);

  // Don't allow to fail
  if (ts.size() == 0) {
    ts.add(Trans3d());
    if (equiv_sets)
      set_equiv_elems_identity(geom, equiv_sets);
  }

  return 1;
}

const char *type_str[Symmetry::Ih + 1] = {
    "Unknown", "C1", "Ci", "Cs", "C",  "Cv", "Ch", "D", "Dv",
    "Dh",      "S",  "T",  "Td", "Th", "O",  "Oh", "I", "Ih"};

Symmetry::Symmetry(const Geometry &geom, vector<vector<set<int>>> *equiv_sets)
{
  init(geom, equiv_sets);
}

Symmetry::Symmetry(int type, int n, const Trans3d &pos, Status *stat)
{
  Status st = init(type, n, pos);
  if (stat)
    *stat = st;
}

Symmetry::Symmetry(const string &name, const Trans3d &pos, Status *stat)
{
  Status st = init(name, pos);
  if (stat)
    *stat = st;
}

Symmetry::Symmetry(const SymmetryAxis &sym_axis, const Vec3d &cent)
{
  init(sym_axis, cent);
}

Status Symmetry::init(const Geometry &geom,
                      vector<vector<set<int>>> *equiv_sets)
{
  sym_type = unknown;
  Transformations ts;
  find_syms(geom, ts, equiv_sets);
  *this = Symmetry(ts);
  return (sym_type != unknown)
             ? Status::ok()
             : Status::error("symmetry type could not be determined");
}

Status Symmetry::init(int type, int n, const Trans3d &pos)
{
  sym_type = unknown;
  nfold = 0;
  axes.clear();
  mirrors.clear();
  sub_syms.clear();
  autos = SymmetryAutos();

  if ((type < C1 || type > Ih))
    return Status::error("unknown symmetry type");
  else if ((type < C || type > S) && n)
    return Status::error("symmetry type doesn't take a number");
  else if ((type >= C && type <= S) && !n)
    return Status::error("symmetry type needs a number");
  else if (type == S && !is_even(n))
    return Status::error("symmetry type S needs an even number");

  to_std = pos;
  sym_type = type;
  if (sym_type >= C && sym_type <= S) // principal axis
    nfold = n;

  // Allow for alternative descriptions of "lesser" symmetries
  if ((nfold < 2 && sym_type >= C && sym_type <= Dh) ||
      (nfold == 2 && sym_type == S)) {
    Symmetry tmp(get_trans());
    *this = tmp;
  }

  return Status::ok();
}

Status Symmetry::init(const string &name, const Trans3d &pos)
{
  sym_type = unknown;
  nfold = 0;

  char numbers[] = "0123456789";
  size_t pos1 = strcspn(name.c_str(), numbers); // pos to end or first digit
  size_t pos2 = pos1;                           // pos to end or first digit
  if (pos1 < name.length())                     // name includes a digit
    pos2 +=
        strspn(name.c_str() + pos1, numbers); // pos to end or after last digit

  // map all the symmetry types to their type number
  map<string, int> sym_names;
  for (int i = 1; i < Ih + 1; i++)
    sym_names[type_str[i]] = i;

  // Assemble the symmetry type string (allows, e.g. Cv3 = C3v, but not 3Cv)
  auto sym_type_str = (pos1) ? name.substr(0, pos1) + name.substr(pos2) : "";
  auto sni = sym_names.find(sym_type_str.c_str());
  int type = (sni != sym_names.end()) ? sni->second : -1;

  int fold = atoi(name.substr(pos1, pos2 - pos1).c_str()); // digit substring

  return init(type, fold, pos);
}

void Symmetry::init(const SymmetryAxis &sym_axis, const Vec3d &cent)
{
  axes.clear();
  mirrors.clear();
  sub_syms.clear();
  autos = SymmetryAutos();

  sym_type = sym_axis.get_sym_type();
  nfold = sym_axis.get_nfold();
  to_std = Trans3d::translate(cent) * Trans3d::align(Vec3d::X, Vec3d::Z,
                                                     sym_axis.get_axis(),
                                                     sym_axis.get_perp());
}

string Symmetry::get_symbol() const
{
  if (sym_type < C || sym_type > S)
    return type_str[sym_type];
  else
    return string() + type_str[sym_type][0] + std::to_string(nfold) +
           (type_str[sym_type] + 1);
}

Transformations &Symmetry::get_trans(Transformations &ts) const
{
  ts.clear();
  if (sym_type >= C1 && sym_type <= Ih) {
    if (sym_type < C || sym_type > S) {
      typedef sch_gen &(sch_gen::*PF_SYM)();
      map<int, PF_SYM> sym_names_other;
      sym_names_other[unknown] = &sch_gen::unit;
      sym_names_other[C1] = &sch_gen::unit;
      sym_names_other[Cs] = &sch_gen::Cs;
      sym_names_other[Ci] = &sch_gen::Ci;
      sym_names_other[T] = &sch_gen::T;
      sym_names_other[Th] = &sch_gen::Th;
      sym_names_other[Td] = &sch_gen::Td;
      sym_names_other[O] = &sch_gen::O;
      sym_names_other[Oh] = &sch_gen::Oh;
      sym_names_other[I] = &sch_gen::I;
      sym_names_other[Ih] = &sch_gen::Ih;

      auto soi = sym_names_other.find(sym_type);
      ts = (sch_gen().*soi->second)();
    }
    else {
      typedef sch_gen &(sch_gen::*PF_SYMN)(int n);
      map<int, PF_SYMN> sym_names_n;
      sym_names_n[C] = &sch_gen::C;
      sym_names_n[Cv] = &sch_gen::Cv;
      sym_names_n[Ch] = &sch_gen::Ch;
      sym_names_n[D] = &sch_gen::D;
      sym_names_n[Dv] = &sch_gen::Dv;
      sym_names_n[Dh] = &sch_gen::Dh;
      sym_names_n[S] = &sch_gen::S;

      auto sni = sym_names_n.find(sym_type);
      ts = (sch_gen().*sni->second)(nfold);
    }
  }

  ts.conjugate(to_std.inverse());
  return ts;
}

Transformations Symmetry::get_trans() const
{
  Transformations ts;
  return get_trans(ts);
}

const set<SymmetryAxis> &Symmetry::get_axes() const
{
  if (axes.size() == 0 && sym_type != unknown) {
    Transformations ts;
    get_trans(ts);
    axes = Symmetry(ts).axes;
  }
  return axes;
}

const set<Vec3d> &Symmetry::get_mirrors() const
{
  if (mirrors.size() == 0 &&
      (sym_type == Cs || sym_type == Cv || sym_type == Ch || sym_type == Dv ||
       sym_type == Dh || sym_type == Td || sym_type == Th || sym_type == Oh ||
       sym_type == Ih)) {
    Transformations ts;
    get_trans(ts);
    for (const auto &t : ts) {
      SymmetryAxis ax = t;
      if (ax.get_sym_type() == Cs)
        mirrors.insert(ax.get_axis());
    }
  }
  return mirrors;
}

bool Symmetry::has_inversion_symmetry() const
{
  if (sym_type == Ih || sym_type == Oh || sym_type == Th || sym_type == Ci ||
      (sym_type == Ch && nfold % 2 == 0) || (sym_type == Dv && nfold % 2) ||
      (sym_type == Dh && nfold % 2 == 0) || (sym_type == S && (nfold / 2) % 2))
    return true;
  else
    return false;
}

void get_equiv_elems(const Geometry &geom, const Transformations &ts,
                     vector<vector<set<int>>> *equiv_sets)
{

  equiv_sets->clear();
  if (ts.size() <= 1) {
    set_equiv_elems_identity(geom, equiv_sets);
    return;
  }

  Geometry merged_geom = geom;
  vector<map<int, set<int>>> orig_equivs;
  merge_coincident_elements(merged_geom, "vef", &orig_equivs, anti::epsilon);
  Geometry test_geom = merged_geom;

  vector<map<int, set<int>>> equiv_elems(3);
  int cnts[3] = {(int)merged_geom.verts().size(),
                 (int)merged_geom.edges().size(),
                 (int)merged_geom.faces().size()};

  for (const auto &t : ts) {
    Geometry trans_geom = merged_geom;
    trans_geom.transform(t);
    vector<map<int, set<int>>> new_equivs;
    check_coincidence(merged_geom, trans_geom, &new_equivs, sym_eps);
    update_equiv_elems(equiv_elems, new_equivs, cnts);
  }

  if (equiv_sets)
    equiv_elems_to_sets(*equiv_sets, equiv_elems, orig_equivs);
}

void Symmetry::add_sub_axes(const Symmetry &sub) const
{
  vector<int> factors;
  int fold = sub.get_nfold();
  for (int i = 1; i <= fold / 2; i++)
    if (fold % i == 0)
      factors.push_back(fold / i);
  for (int factor : factors) {
    int nfold = factor;
    Symmetry sub_ax(sub.get_sym_type(), nfold, sub.get_to_std());
    switch (sub.get_sym_type()) {
    case C:
      sub_syms.insert(sub_ax);
      break;
    case Cv:
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(C);
      sub_syms.insert(sub_ax);
      if (fold % 2 == 0) {
        sub_ax.set_sym_type(Cv);
        sub_ax.set_to_std(Trans3d::rotate(0, 0, M_PI / (fold)) *
                          sub.get_to_std());
        sub_syms.insert(sub_ax);
      }
      break;
    case Ch:
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(C);
      sub_syms.insert(sub_ax);
      if (factor > 2 && factor % 2 == 0) {
        sub_ax.set_sym_type(S);
        sub_syms.insert(sub_ax);
      }
      break;
    case D:
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(C);
      sub_syms.insert(sub_ax);
      break;
    case Dv:
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(D);
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(C);
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(Cv);
      sub_ax.set_to_std(Trans3d::rotate(0, 0, M_PI / (2 * nfold)) *
                        sub.get_to_std());
      sub_syms.insert(sub_ax);
      if (nfold % 2 == 0) {
        sub_ax.set_to_std(Trans3d::rotate(0, 0, -M_PI / (2 * nfold)) *
                          sub.get_to_std());
        sub_syms.insert(sub_ax);
      }
      sub_syms.insert(sub_ax);
      sub_ax.set_nfold(2 * nfold);
      sub_ax.set_to_std(sub.get_to_std());
      sub_ax.set_sym_type(S);
      sub_syms.insert(sub_ax);
      break;
    case Dh:
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(D);
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(C);
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(Cv);
      sub_syms.insert(sub_ax);
      sub_ax.set_sym_type(Ch);
      sub_syms.insert(sub_ax);
      if (nfold > 2 && nfold % 2 == 0) {
        sub_ax.set_sym_type(S);
        sub_syms.insert(sub_ax);
        sub_ax.set_nfold(nfold / 2);
        sub_ax.set_sym_type(Dv);
        sub_syms.insert(sub_ax);
      }
      break;
    case S:
      if (nfold > 2 && nfold % 2 == 0) {
        sub_syms.insert(sub_ax);
        sub_ax.set_sym_type(C);
        sub_ax.set_nfold(nfold / 2);
        sub_syms.insert(sub_ax);
      }
      break;
    }
  }
}

class same_sym_group {
private:
  const Symmetry &sym;

public:
  same_sym_group(const Symmetry &s) : sym(s) {}

  bool operator()(const Symmetry &cmp) const
  {
    if (cmp.get_sym_type() == sym.get_sym_type()) {
      if (sym.get_sym_type() < Symmetry::C || sym.get_sym_type() > Symmetry::S)
        return true;
      else if (cmp.get_nfold() == sym.get_nfold())
        return true;
    }
    return false;
  }
};

static Status create_sub_sym(const Symmetry &full, Symmetry *sub,
                             const Symmetry &sub_sym, int conj_type)
{
  *sub = Symmetry(Symmetry::unknown);
  if (conj_type < 0)
    return Status::error("conjugation type number cannot be negative");

  const set<Symmetry> &sub_syms = full.get_sub_syms();
  auto si = sub_syms.begin();
  same_sym_group cmp(sub_sym);

  // To minimise realignment of first conjugate subgroup, order on magnitude
  // of difference of to_std rotations applied to a general test vector
  std::multimap<double, Symmetry> conj_syms;
  while ((si = find_if(si, sub_syms.end(), cmp)) != sub_syms.end()) {
    const Vec3d test_v(1.1, 2.0, M_PI); // not on any axis!
    Vec3d d = full.get_to_std() * test_v - si->get_to_std() * test_v;
    conj_syms.insert(std::make_pair(d.len(), *si));
    ++si;
  }

  if (!conj_syms.size())
    return Status::error(msg_str("%s is not a sub-symmetry of %s",
                                 sub_sym.get_symbol().c_str(),
                                 full.get_symbol().c_str()));
  if (conj_type >= (int)conj_syms.size())
    return Status::error(
        msg_str("conjugation type too large for %s (last number: %d)",
                sub_sym.get_symbol().c_str(), (int)conj_syms.size() - 1));

  // valid
  auto mi = conj_syms.begin();
  int cnt = 0;
  while (cnt++ != conj_type)
    mi++;
  *sub = mi->second;

  return Status::ok();
}

Status Symmetry::get_sub_sym(const Symmetry &sub_sym, Symmetry *sub,
                             int conj_type) const
{
  return create_sub_sym(*this, sub, sub_sym, conj_type);
}

Status Symmetry::get_sub_sym(const string &sub_name, Symmetry *sub) const
{
  *sub = Symmetry();

  Split parts(sub_name, ",");
  if (parts.size() > 2)
    return Status::error("too many comma separated parts");

  Status stat;
  Symmetry sub_sym;
  if (parts.size() == 0 || strncmp(parts[0], "full", strlen(parts[0])) == 0)
    sub_sym = *this;
  else if (!(stat = sub_sym.init(parts[0], Trans3d())))
    return Status::error(msg_str("sub-symmetry type: %s", stat.c_msg()));

  int sub_sym_conj = 0;
  if (parts.size() > 1 && !(stat = read_int(parts[1], &sub_sym_conj)))
    return Status::error(
        msg_str("sub-symmetry conjugation number: %s", stat.c_msg()));

  if (!(stat = get_sub_sym(sub_sym, sub, sub_sym_conj)))
    return Status::error(msg_str("sub-symmetry: %s", stat.c_msg()));

  return Status::ok();
}

const set<Symmetry> &Symmetry::get_sub_syms() const
{
  sub_syms.clear();
  const Vec3d axis = Vec3d::Z;
  const Vec3d perp = Vec3d::X;
  int dih_type;
  if (sub_syms.size() == 0 && sym_type != unknown) {
    Symmetry sym = *this;
    sub_syms.insert(sym);
    if (sym_type != C1) { // Don't add C1 a second time
      sym.sym_type = C1;
      sub_syms.insert(sym);
    }
    if (has_inversion_symmetry()) {
      sym.sym_type = Ci;
      sub_syms.insert(sym);
    }
    switch (sym_type) {
    case Ih:
      sym.sym_type = I;
      sub_syms.insert(sym);
      sym.sym_type = Th;
      sub_syms.insert(sym);
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.sym_type = Ci;
      sub_syms.insert(sym);
      sym.sym_type = Cs;
      sub_syms.insert(sym);
      sym.init(Dv, 5,
               Trans3d::align(A5, Vec3d::Z, axis,
                              Trans3d::rotate(axis, M_PI / 10) * Vec3d::X) *
                   to_std);
      add_sub_axes(sym);
      sym.init(Dv, 3,
               Trans3d::align(A3, A5, axis,
                              Trans3d::rotate(axis, M_PI / 6) * Vec3d::X) *
                   to_std);
      add_sub_axes(sym);
      sym.init(Dh, 2, to_std);
      add_sub_axes(sym);
      break;

    case I:
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.init(D, 5,
               Trans3d::align(A5, Vec3d::Z, axis,
                              Trans3d::rotate(axis, M_PI / 10) * Vec3d::X) *
                   to_std);
      add_sub_axes(sym);
      sym.init(D, 3,
               Trans3d::align(A3, A5, axis,
                              Trans3d::rotate(axis, M_PI / 6) * Vec3d::X) *
                   to_std);
      add_sub_axes(sym);
      sym.init(D, 2, to_std);
      add_sub_axes(sym);
      break;

    case Oh:
      sym.sym_type = O;
      sub_syms.insert(sym);
      sym.sym_type = Td;
      sub_syms.insert(sym);
      sym.sym_type = Th;
      sub_syms.insert(sym);
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.sym_type = Ci;
      sub_syms.insert(sym);
      sym.sym_type = Cs;
      sub_syms.insert(sym);
      sym.init(Cs, 0, Trans3d::rotate(Vec3d(0, 1, 1), axis) * to_std);
      add_sub_axes(sym);
      sym.init(Dh, 4, to_std);
      add_sub_axes(sym);
      sym.init(Dv, 3, Trans3d::align(A3, Vec3d(1, -1, 0), axis, perp) * to_std);
      add_sub_axes(sym);
      sym.init(Dh, 2, Trans3d::rotate(Vec3d(0, 1, 1), axis) * to_std);
      add_sub_axes(sym);
      break;

    case O:
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.init(D, 4, to_std);
      add_sub_axes(sym);
      sym.init(D, 3, Trans3d::align(A3, Vec3d(1, -1, 0), axis, perp) * to_std);
      add_sub_axes(sym);
      sym.init(D, 2, Trans3d::rotate(Vec3d(0, 1, 1), axis) * to_std);
      add_sub_axes(sym);
      break;

    case Th:
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.sym_type = Ci;
      sub_syms.insert(sym);
      sym.sym_type = Cs;
      sub_syms.insert(sym);
      sym.init(S, 6, Trans3d::rotate(A3, axis) * to_std);
      add_sub_axes(sym);
      sym.init(Dh, 2, to_std);
      add_sub_axes(sym);
      break;

    case Td:
      sym.sym_type = T;
      sub_syms.insert(sym);
      sym.init(Cs, 0, Trans3d::rotate(Vec3d(0, 1, 1), axis) * to_std);
      sub_syms.insert(sym);
      sym.init(Cv, 3,
               Trans3d::align(A3, Vec3d(1, -1, -1), axis, perp) * to_std);
      add_sub_axes(sym);
      sym.init(Dv, 2, to_std);
      add_sub_axes(sym);
      break;

    case T:
      sym.init(C, 3, Trans3d::rotate(A3, axis) * to_std);
      add_sub_axes(sym);
      sym.init(D, 2, to_std);
      add_sub_axes(sym);
      break;

    case Dh:
      add_sub_axes(*this);
      sym.init(Cs, 0, to_std); // horizontal mirror
      sub_syms.insert(sym);
      if (nfold % 2 == 0) { // nfold even: second axis and vert mirror
        dih_type = Dh;
        // do this first to prefer main axis D groups
        if (nfold > 2) {
          sym.init(Dh, nfold / 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
          add_sub_axes(sym);
          sym.init(Dv, nfold / 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
          add_sub_axes(sym);
        }
        // vertical mirror through dihedral axis
        sym.init(Cs, 0,
                 Trans3d::rotate(Vec3d::X, M_PI / 2) *
                     Trans3d::rotate(axis, M_PI / nfold) * to_std);
        sub_syms.insert(sym);
        // dihedral axis
        sym.init(dih_type, 2,
                 Trans3d::rotate(Vec3d::Y, M_PI / 2) *
                     Trans3d::rotate(axis, M_PI / nfold) * to_std);
        add_sub_axes(sym);
      }
      else
        dih_type = Cv;

      // vertical mirror through dihedral axis
      sym.init(Cs, 0, Trans3d::rotate(Vec3d::X, M_PI / 2) * to_std);
      sub_syms.insert(sym);

      sym.init(dih_type, 2, Trans3d::rotate(Vec3d::Y, M_PI / 2) * to_std);
      add_sub_axes(sym);
      break;

    case Dv:
      add_sub_axes(*this);
      // vertical mirror between dihedral axes
      sym.init(Cs, 0,
               Trans3d::rotate(axis, M_PI / (2 * nfold)) *
                   Trans3d::rotate(Vec3d::X, M_PI / 2) * to_std);
      sub_syms.insert(sym);
      if (nfold % 2 == 0) {
        dih_type = D;
        if (nfold > 2) {
          sym.init(Dv, nfold / 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
          add_sub_axes(sym);
          // do this first to prefer main axis D groups
          sym.init(D, 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
          add_sub_axes(sym);
        }
        sym.init(dih_type, 2,
                 Trans3d::rotate(Vec3d::Y, M_PI / 2) *
                     Trans3d::rotate(axis, M_PI / nfold) * to_std);
        add_sub_axes(sym);
      }
      else
        dih_type = Ch;
      sym.init(dih_type, 2, Trans3d::rotate(Vec3d::Y, M_PI / 2) * to_std);
      add_sub_axes(sym);
      break;

    case D:
      add_sub_axes(*this);
      if (nfold % 2 == 0 && nfold > 2) {
        sym.init(D, nfold / 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
        add_sub_axes(sym);
      }
      dih_type = (nfold % 2 == 0) ? D : C;

      sym.init(dih_type, 2,
               Trans3d::rotate(Vec3d::Y, M_PI / 2) *
                   Trans3d::rotate(axis, M_PI / nfold) * to_std);
      add_sub_axes(sym);
      sym.init(dih_type, 2, Trans3d::rotate(Vec3d::Y, M_PI / 2) * to_std);
      add_sub_axes(sym);
      break;

    case S:
      add_sub_axes(*this);
      break;

    case Ch:
      sym.sym_type = Cs; // horizontal mirror
      sub_syms.insert(sym);
      add_sub_axes(*this);
      break;

    case Cv:
      if (nfold % 2 == 0) { // nfold even: add second vertical mirror
        sym.init(Cs, 0,
                 Trans3d::rotate(axis, M_PI / nfold) *
                     Trans3d::rotate(Vec3d::X, M_PI / 2) * to_std);
        sub_syms.insert(sym);
        if (nfold > 2) {
          sym.init(Cv, nfold / 2, Trans3d::rotate(axis, M_PI / nfold) * to_std);
          sub_syms.insert(sym);
        }
      }
      // vertical mirror
      sym.init(Cs, 0, Trans3d::rotate(Vec3d::X, M_PI / 2) * to_std);
      sub_syms.insert(sym);
      add_sub_axes(*this);
      break;

    case C:
      add_sub_axes(*this);
      break;

    // No proper subgroups
    case Ci:
    case Cs:
    case C1:
      break;
    }
  }

  return sub_syms;
}

Symmetry Symmetry::get_max_direct_sub_sym() const
{
  int type;
  int fold = 0;
  if (sym_type == Ih || sym_type == I)
    type = I;
  else if (sym_type == Oh || sym_type == O)
    type = O;
  else if (sym_type == Th || sym_type == Td || sym_type == T)
    type = T;
  else if (sym_type == Dh || sym_type == Dv || sym_type == D) {
    type = D;
    fold = nfold;
  }
  else if (sym_type == S || sym_type == Ch || sym_type == Cv || sym_type == C) {
    type = C;
    fold = (sym_type == S) ? nfold / 2 : nfold;
  }
  else
    type = C1;

  return Symmetry(type, fold, to_std);
}

SymmetryAutos &Symmetry::get_autos()
{
  if (!autos.is_set())
    autos = SymmetryAutos(*this);
  return autos;
}

Subspace Symmetry::get_fixed_subspace() const
{
  Subspace::SubspaceType type;
  Vec3d point = get_to_std().inverse() * Vec3d(0, 0, 0);
  Vec3d direction;

  switch (get_sym_type()) {
  case Symmetry::C1:
    type = Subspace::SubspaceType::space;
    break;
  case Symmetry::Cs:
    type = Subspace::SubspaceType::plane;
    if (get_mirrors().begin() != get_mirrors().end())
      direction = *get_mirrors().begin(); // mirror normal
    break;
  case Symmetry::C:
  case Symmetry::Cv:
    type = Subspace::SubspaceType::line;
    if (get_axes().begin() != get_axes().end())
      direction = get_axes().begin()->get_axis();
    break;
  case Symmetry::Ci:
  case Symmetry::S:
  case Symmetry::Ch:
  case Symmetry::D:
  case Symmetry::Dv:
  case Symmetry::Dh:
  case Symmetry::T:
  case Symmetry::Td:
  case Symmetry::Th:
  case Symmetry::O:
  case Symmetry::Oh:
  case Symmetry::I:
  case Symmetry::Ih:
    type = Subspace::SubspaceType::point;
    break;
  default:
    type = Subspace::SubspaceType::none;
  }
  return Subspace(type, point, direction);
}

bool Symmetry::operator<(const Symmetry &s) const
{
  if (sym_type < s.sym_type)
    return true;
  else if (sym_type > s.sym_type)
    return false;

  if (sym_type >= C && sym_type <= S) {
    if (nfold < s.nfold)
      return true;
    else if (nfold > s.nfold)
      return false;
  }

  // Expensive test if same kind of symmetry group!!!
  return compare(get_trans(), s.get_trans()) == -1;
}

void SymmetryAutos::init()
{
  free_vars = FREE_NONE;
  fixed_type = 0;
  fixed_trans.clear();
  for (int i = 0; i < 3; i++)
    rot[i] = transl[i] = 0.0;
}

SymmetryAutos::SymmetryAutos(const Symmetry &sym)
{
  init();
  int type = sym.get_sym_type();
  if (type == Symmetry::unknown) // fixed_trans empty, indicates unset;
    return;

  int nfold = sym.get_nfold();

  vector<vector<Trans3d>> fixed;       // fixed (non-free) realignments
  fixed.push_back(vector<Trans3d>(1)); // identity

  if (type == Symmetry::I || type == Symmetry::O || type == Symmetry::T ||
      type == Symmetry::D || type == Symmetry::C || type == Symmetry::C1) {
    fixed.push_back(vector<Trans3d>(2));
    fixed.back()[1] = Trans3d::inversion(); // reflect through origin
  }

  if (type == Symmetry::Td || type == Symmetry::T || type == Symmetry::S ||
      type == Symmetry::Ch || type == Symmetry::Cv || type == Symmetry::C) {
    fixed.push_back(vector<Trans3d>(2));
    fixed.back()[1] = Trans3d::rotate(M_PI, 0, 0); // flip "principal" axis
  }

  if (type == Symmetry::Dh || type == Symmetry::Dv || type == Symmetry::D ||
      type == Symmetry::Cv) {
    fixed.push_back(vector<Trans3d>(2));
    fixed.back()[1] =
        Trans3d::rotate(0, 0, M_PI / nfold); // rotate base vert to edge
  }

  if ((type == Symmetry::Dh || type == Symmetry::D) && nfold == 2) {
    fixed.push_back(vector<Trans3d>(3));
    fixed.back()[1] =
        Trans3d::rotate(Vec3d(1, 1, 1), 2 * M_PI / 3); // rotate D2 axes
    fixed.back()[2] =
        Trans3d::rotate(Vec3d(1, 1, 1), -2 * M_PI / 3); // rotate D2 axes
  }

  // find all combinations of transformations involving one member
  // from each set (is there an STL algorithm for this?)
  int sz = fixed.size();
  vector<unsigned int> visit_idxs(sz, 0);
  while (visit_idxs[0] < fixed[0].size()) { // until first counter rolls
    Trans3d trans;
    for (int i = 0; i < sz; i++)
      trans = fixed[i][visit_idxs[i]] * trans;
    fixed_trans.push_back(trans);
    // Trans3d::inverse(sym.get_to_std()) * trans * sym.get_to_std() );

    // increment visit counters
    visit_idxs[sz - 1]++;
    for (int j = 1; j < sz; j++)
      if (visit_idxs[sz - j] >= fixed[sz - j].size()) {
        visit_idxs[sz - j] = 0;
        visit_idxs[sz - j - 1]++;
      }
  }

  if (type == Symmetry::S || type == Symmetry::Ch)
    free_vars = FREE_ROT_PRINCIPAL;
  if (type == Symmetry::Cv)
    free_vars = FREE_TRANSL_PRINCIPAL;
  else if (type == Symmetry::C)
    free_vars = FREE_ROT_PRINCIPAL | FREE_TRANSL_PRINCIPAL;
  else if (type == Symmetry::Cs)
    free_vars = FREE_ROT_PRINCIPAL | FREE_TRANSL_PLANE;
  else if (type == Symmetry::Ci)
    free_vars = FREE_ROT_FULL;
  else if (type == Symmetry::C1)
    free_vars = FREE_ROT_FULL | FREE_TRANSL_SPACE;
}

void SymmetryAutos::set_fixed(const Transformations &fixed)
{
  fixed_trans.clear();
  vector<Trans3d> indirect;
  fixed_trans.push_back(Trans3d());
  set<Trans3d>::const_iterator si;
  for (si = fixed.begin(); si != fixed.end(); ++si) {
    Isometry iso(*si);
    if (iso.is_direct() && iso.get_rot_type() != Isometry::rt_unit)
      fixed_trans.push_back(*si);
    else if (!iso.is_direct())
      indirect.push_back(*si);
  }
  fixed_trans.insert(fixed_trans.end(), indirect.begin(), indirect.end());
}

Status SymmetryAutos::set_fixed_type(int type)
{
  if (!is_set())
    return Status::error("fixed type: symmetry type not yet set");

  if (type < 0 || type >= (int)fixed_trans.size()) {
    string range = (fixed_trans.size() == 1)
                       ? "0"
                       : msg_str("0-%d", (int)fixed_trans.size() - 1);
    return Status::error(msg_str("fixed type: invalid type %d, must be %s",
                                 type, range.c_str()));
  }

  fixed_type = type;

  return Status::ok();
}

Status SymmetryAutos::set_rot_principal(double rot_ang)
{
  if (!(free_vars & FREE_ROT_PRINCIPAL))
    return Status::error(
        "symmetry type does not have a free principal axis rotation");

  rot[0] = deg2rad(rot_ang);
  rot[1] = 0.0;
  rot[2] = 0.0;

  return Status::ok();
}

Status SymmetryAutos::set_rot_full(double rot_x, double rot_y, double rot_z)
{
  if (!(free_vars & FREE_ROT_FULL))
    return Status::error("symmetry type does not have a free full rotation");

  rot[0] = deg2rad(rot_x);
  rot[1] = deg2rad(rot_y);
  rot[2] = deg2rad(rot_z);

  return Status::ok();
}

Status SymmetryAutos::set_transl_principal(double transl0)
{
  if (!(free_vars & FREE_TRANSL_PRINCIPAL))
    return Status::error(
        "symmetry type does not have a free principal axis translation");

  transl[0] = transl0;
  transl[1] = 0.0;
  transl[2] = 0.0;

  return Status::ok();
}

Status SymmetryAutos::set_transl_plane(double transl0, double transl1)
{
  if (!(free_vars & FREE_TRANSL_PLANE))
    return Status::error(
        "symmetry type does not have a free plane translation");

  transl[0] = transl0;
  transl[1] = transl1;
  transl[2] = 0.0;

  return Status::ok();
}

Status SymmetryAutos::set_transl_space(double transl0, double transl1,
                                       double transl2)
{
  if (!(free_vars & FREE_TRANSL_SPACE))
    return Status::error(
        "symmetry type does not have a free space translation");

  transl[0] = transl0;
  transl[1] = transl1;
  transl[2] = transl2;

  return Status::ok();
}

int SymmetryAutos::num_free_rots() const
{
  int cnt = 0;
  if (free_vars & FREE_ROT_PRINCIPAL)
    cnt = 1;
  else if (free_vars & FREE_ROT_FULL)
    cnt = 3;
  return cnt;
}

int SymmetryAutos::num_free_transls() const
{
  int cnt = 0;
  if (free_vars & FREE_TRANSL_PRINCIPAL)
    cnt = 1;
  else if (free_vars & FREE_TRANSL_PLANE)
    cnt = 2;
  else if (free_vars & FREE_TRANSL_SPACE)
    cnt = 3;
  return cnt;
}

Trans3d SymmetryAutos::get_realignment() const
{
  Trans3d trans;
  if (is_set())
    trans = fixed_trans[fixed_type];

  if (free_vars & FREE_ROT_PRINCIPAL)
    trans = Trans3d::rotate(Vec3d::Z, rot[0]) * trans;
  else if (free_vars & FREE_ROT_FULL)
    trans = Trans3d::rotate(rot[0], rot[1], rot[2]) * trans;

  if (free_vars & FREE_TRANSL_PRINCIPAL) // z-axis
    trans = Trans3d::translate(Vec3d(0, 0, transl[0])) * trans;
  else if (free_vars & FREE_TRANSL_PLANE) // xy-plane
    trans = Trans3d::translate(Vec3d(transl[0], transl[1], 0)) * trans;
  else if (free_vars & FREE_TRANSL_SPACE)
    trans = Trans3d::translate(Vec3d(transl[0], transl[1], transl[2])) * trans;

  return trans;
}

Status SymmetryAutos::set_realignment(const string &realign)
{
  vector<double> vars;

  string realign_cpy(realign);         // copy, do not access as C++ string
  char *realign_ptr = &realign_cpy[0]; // may be used to modify characters

  char *p = strchr(realign_ptr, ':');
  if (p)
    *p = '\0'; // terminator at first ':'

  Status stat;
  int type = 0;
  if (realign == "" || (p && p == realign_ptr)) // empty, so set default of 0
    type = 0;
  else if (!(stat = read_int(realign_ptr, &type)))
    return Status::error(msg_str("fixed type: %s", stat.c_msg()));

  if (!(stat = set_fixed_type(type)))
    return stat;

  if (p && *(p + 1)) { // colon found and characters afterwards
    if (!(stat = read_double_list(p + 1, vars, 0, ":")))
      return Status::error(msg_str("free variables: %s", stat.c_msg()));
  }

  int rot_cnt = num_free_rots();       // number of rotation angles
  int transl_cnt = num_free_transls(); // number of translation distances

  int vars_sz = vars.size();
  if (vars_sz && vars_sz != rot_cnt && vars_sz != rot_cnt + transl_cnt) {
    string msg = msg_str("free variable list: %d numbers given", vars.size());
    if (rot_cnt + transl_cnt == 0)
      msg += " but there are no free variables";
    else if (rot_cnt > 0 && transl_cnt > 0)
      msg += msg_str(", must give %d (rotation) or %d (rotation then "
                     "translation) colon separated numbers",
                     rot_cnt, rot_cnt + transl_cnt);
    else if (rot_cnt > 0)
      msg +=
          msg_str(", must give %d (rotation) colon separated numbers", rot_cnt);
    else //(transl_cnt>0)
      msg += msg_str(", must give %d (translation) colon separated numbers",
                     transl_cnt);

    return Status::error(msg);
  }

  if (vars_sz) {
    if (free_vars & FREE_ROT_PRINCIPAL)
      set_rot_principal(vars[0]);
    else if (free_vars & FREE_ROT_FULL)
      set_rot_full(vars[0], vars[1], vars[2]);
  }

  if (vars_sz > rot_cnt) {
    if (free_vars & FREE_TRANSL_PRINCIPAL)
      set_transl_principal(vars[rot_cnt + 0]);
    else if (free_vars & FREE_TRANSL_PLANE)
      set_transl_plane(vars[rot_cnt + 0], vars[rot_cnt + 1]);
    else if (free_vars & FREE_TRANSL_SPACE)
      set_transl_space(vars[rot_cnt + 0], vars[rot_cnt + 1], vars[rot_cnt + 2]);
  }

  return Status::ok();
}

//-----------------------------------------------------------------------
// Orbits and Stabilzers

static void copy_face_to_geom(const Geometry &geom, int f_idx, Geometry *o_geom)
{
  o_geom->clear_all();
  const int f_sz = geom.faces(f_idx).size();
  vector<int> f_new(f_sz);
  map<int, int> v_orig2v_new;
  for (int i = 0; i < f_sz; i++) {
    const int v_orig = geom.faces(f_idx, i); // original index number
    int v_new;                               // new vertex number
    map<int, int>::iterator mi = v_orig2v_new.find(geom.faces(f_idx, i));
    if (mi == v_orig2v_new.end()) {
      v_new = o_geom->verts().size();      // new map, use next integer
      v_orig2v_new[v_orig] = v_new;        // add map
      o_geom->add_vert(geom.verts(v_new)); // copy vertex, ignore colours
    }
    else
      v_new = mi->second; // index previously mapped
    f_new[i] = v_new;     // use mapped vertex in face
  }
  o_geom->add_face(f_new); // (ignore colours)
}

static Symmetry get_elem_stabilizer(const Geometry &geom,
                                    const vector<int> &elem,
                                    const Symmetry &sym)
{
  if (elem.size() == 0)
    return Symmetry();

  const Transformations &ts = sym.get_trans();
  Transformations stab_trans;
  const int elem_sz = elem.size();
  vector<Vec3d> trans_elem(elem_sz); // transformed element as cycle of coords
  for (set<Trans3d>::iterator si = ts.get_trans().begin();
       si != ts.get_trans().end(); si++) {
    for (int i = 0; i < elem_sz; i++)
      trans_elem[i] = (*si) * geom.verts(elem[i]);

    // starting at each vertex, see if cycle is fully coincident
    bool coincident = false;
    for (int i = 0; i < elem_sz; i++) { // index
      coincident = true;
      for (int off = 0; off < elem_sz; off++) { // offset
        if (compare(geom.verts(elem[i]), trans_elem[(i + off) % elem_sz],
                    anti::epsilon) != 0) {
          coincident = false;
          break;
        }
      }
      if (coincident)
        break;
    }

    if (coincident) // element carried onto itself
      stab_trans += *si;
  }

  Symmetry stab(stab_trans);
  return stab;
}

static Symmetry get_face_stabilizer(const Geometry &geom, int f_idx,
                                    const Symmetry &sym)
{
  Geometry f_geom, test_geom;
  copy_face_to_geom(geom, f_idx, &f_geom);

  Transformations stab_trans;
  const Transformations &ts = sym.get_trans();
  for (set<Trans3d>::iterator si = ts.get_trans().begin();
       si != ts.get_trans().end(); si++) {
    test_geom = f_geom;
    test_geom.transform(*si);
    test_geom.append(f_geom);
    merge_coincident_elements(test_geom, "f", anti::epsilon);
    if (test_geom.faces().size() == 1) // face carried onto itself
      stab_trans += *si;
  }

  Symmetry stab(stab_trans);
  return stab;
}

static Symmetry get_edge_stabilizer(const Geometry &geom, int e_idx,
                                    const Symmetry &sym)
{
  return get_elem_stabilizer(geom, geom.edges(e_idx), sym);
}

static Symmetry get_vert_stabilizer(const Geometry &geom, int v_idx,
                                    const Symmetry &sym)
{
  Transformations stab_trans;
  const Transformations &ts = sym.get_trans();
  for (set<Trans3d>::iterator si = ts.get_trans().begin();
       si != ts.get_trans().end(); si++) {
    if (compare(geom.verts(v_idx), (*si) * geom.verts(v_idx), sym_eps) == 0)
      stab_trans += *si;
  }

  Symmetry stab(stab_trans);
  return stab;
}

Symmetry get_elem_stabilizer(const Geometry &geom, int idx, int elem_type,
                             const Symmetry &sym)
{
  if (elem_type == ELEM_VERTS)
    return get_vert_stabilizer(geom, idx, sym);
  else if (elem_type == ELEM_EDGES)
    return get_edge_stabilizer(geom, idx, sym);
  else if (elem_type == ELEM_FACES)
    return get_face_stabilizer(geom, idx, sym);
  else // called with invalid type
    return Symmetry();
}

bool elem_fixed(const Geometry &geom, int idx, int elem_type,
                const Symmetry &sym)
{
  Symmetry stab = get_elem_stabilizer(geom, idx, elem_type, sym);
  return stab.get_trans() == sym.get_trans();
}

//-----------------------------------------------------------------
// SymmetricUpdater

// Initialise vertex orbit
void SymmetricUpdater::init_vert_orbit(int orbit_idx, const set<int> &orbit)
{
  // fprintf(stderr, "\ninit_vert_orbit in\n");
  vector<int> idxs(orbit.begin(), orbit.end());
  int v_idx = idxs[0];
  // fprintf(stderr, "v_idx = %d\n", v_idx);
  Symmetry stab = get_vert_stabilizer(geom, v_idx, symmetry);
  // fprintf(stderr, "stab=%s: ", stab.get_symbol().c_str());

  // Set of transformations that carry the vertex once onto each orbit vertex
  Transformations trans;
  trans.min_set(symmetry.get_trans(), stab.get_trans());
  // fprintf(stderr, "trans.size()=%lu: ", (unsigned long)trans.size());
  const int num_repeats = trans.size();

  Geometry merge_geom, part_geom;
  part_geom.add_vert(geom.verts(v_idx));
  sym_repeat(merge_geom, part_geom, trans);
  for (auto idx : idxs)
    merge_geom.add_vert(geom.verts(idx));

  vector<map<int, set<int>>> equivs;
  merge_coincident_elements(merge_geom, "v", &equivs, sym_eps);

  // map first of any coincident vertices to its trans order number
  vector<int> idxs_in_trans_order(num_repeats);

  // Debug print
  if (false) {
    for (auto mi = equivs[VERTS].begin(); mi != equivs[VERTS].end(); ++mi) {
      fprintf(stderr, "idx=%d: ", mi->first);
      for (auto si = mi->second.begin(); si != mi->second.end(); si++) {
        if (*si < (int)idxs.size())
          fprintf(stderr, "%d  ", idxs[*si]);
        else
          fprintf(stderr, "*%d", *si);
      }
      fprintf(stderr, "\n");
    }
  }

  for (map<int, set<int>>::iterator mi = equivs[VERTS].begin();
       mi != equivs[VERTS].end(); ++mi) {
    if (mi->second.size() < 2) { // Shouldn't happen!
      ; // fprintf(stderr, "orbit: missing map! (for idx=%d)\n", mi->first);
      // return;
    }
    else {
      set<int>::iterator si = mi->second.begin();
      const int pos = *si;
      const int idx = *(++si);
      idxs_in_trans_order[pos] = idxs[idx - num_repeats];
    }
  }

  map<int, Trans3d> elems;
  int pos = 0;
  for (auto si = trans.begin(); si != trans.end(); si++) {
    // fprintf(stderr, "pos %d: idx %d\n", pos, idxs_in_trans_order[pos]);
    // si->dump();
    elems[idxs_in_trans_order[pos]] = *si;
    pos++;
  }

  // Update orbit information
  orbit_vertex_idx.push_back(v_idx);
  orbit_invariant_subspaces.push_back(stab.get_fixed_subspace());
  for (const auto &elem : elems) {
    // fprintf(stderr, "elem: idx=%d\n", elem.first);
    const auto it_to = transformations.get_trans().find(elem.second);
    const auto it_from =
        transformations.get_trans().find(elem.second.inverse());
    orbit_mapping[elem.first] = ElemOrbitMapping(orbit_idx, it_from, it_to);
  }
  // fprintf(stderr, "init_vert_orbit out\n");
}

SymmetricUpdater::SymmetricUpdater(const Geometry &base_geom, Symmetry sym)
    : symmetry(sym)
{
  geom = base_geom;
  orbit_vertex_idx.clear();
  orbit_invariant_subspaces.clear();
  orbit_mapping.clear();

  transformations = symmetry.get_trans();
  orbit_mapping.resize(
      geom.verts().size(),
      ElemOrbitMapping(-1, transformations.end(), transformations.end()));
  get_equiv_elems(geom, transformations, &equiv_sets);
  const vector<set<int>> &v_equiv_sets = equiv_sets[VERTS];
  for (size_t orbit_idx = 0; orbit_idx < v_equiv_sets.size(); orbit_idx++)
    init_vert_orbit(orbit_idx, v_equiv_sets[orbit_idx]);
}

vector<int> SymmetricUpdater::get_principal(int type)
{
  vector<int> principal_idxs;
  if (type >= VERTS && type <= FACES) {
    for (auto &orbit : get_equiv_sets(type))
      principal_idxs.push_back(*orbit.begin());
  }

  return principal_idxs;
}

void SymmetricUpdater::to_unique_index_list(vector<int> &idxs)
{
  sort(idxs.begin(), idxs.end());
  idxs.erase(unique(idxs.begin(), idxs.end()), idxs.end());
}

vector<int> SymmetricUpdater::sequential_index_list(int size)
{
  vector<int> idxs;
  if (size > 0) {
    idxs.resize(size);
    std::iota(idxs.begin(), idxs.end(), 0);
  }
  return idxs;
}

vector<int>
SymmetricUpdater::get_included_verts(const vector<int> &elem_idxs,
                                     const vector<vector<int>> &elem_verts)
{
  vector<int> included_verts;
  for (int elem_idx : elem_idxs) {
    included_verts.insert(included_verts.end(), elem_verts[elem_idx].begin(),
                          elem_verts[elem_idx].end());
  }
  to_unique_index_list(included_verts);
  return included_verts;
}

vector<int>
SymmetricUpdater::get_associated_elems(const vector<vector<int>> &elems)
{
  vector<int> elem_idxs;
  for (auto &v_orbit : get_equiv_sets(VERTS)) {
    int v_idx = *v_orbit.begin();
    elem_idxs.insert(elem_idxs.end(), elems[v_idx].begin(), elems[v_idx].end());
  }
  to_unique_index_list(elem_idxs);
  return elem_idxs;
}

void SymmetricUpdater::update_principal_vertex(int v_idx, Vec3d point)
{
  const auto &orb = orbit_mapping[v_idx];
  int orb_vert_idx = orbit_vertex_idx[orb.get_orbit_no()];

  // Update the principal orbit vertex with the new value
  if (orb_vert_idx != v_idx)
    geom.verts(orb_vert_idx) = orb.get_trans_to() * point;
  else // don't map if this is the principal vertex
    geom.verts(orb_vert_idx) = point;
  const auto &subspace = orbit_invariant_subspaces[orb.get_orbit_no()];
  geom.verts(orb_vert_idx) =
      subspace.nearest_point(geom.verts(orb_vert_idx)).with_len(point.len());
}

Vec3d SymmetricUpdater::update_from_principal_vertex(int v_idx)
{
  const auto &orb = orbit_mapping[v_idx];
  const int orb_vert_idx = orbit_vertex_idx[orb.get_orbit_no()];

  // Update the value from the principal orbit vertex
  if (orb_vert_idx != v_idx) // don't map if this is the principal vertex
    geom.verts(v_idx) = orb.get_trans_from() * geom.verts(orb_vert_idx);

  return geom.verts(v_idx);
}

void SymmetricUpdater::update_all()
{
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    update_from_principal_vertex(i);
}

const Geometry &SymmetricUpdater::get_geom_final()
{
  update_all();
  return geom;
}

} // namespace anti
