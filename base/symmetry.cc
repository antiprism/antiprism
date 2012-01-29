/*
   Copyright (c) 2008, Adrian Rossiter

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


#include <string.h>
#include <stdlib.h>
#include <set>
#include <algorithm>

#include "utils.h"
#include "math_utils.h"
#include "transforms.h"
#include "info.h"
#include "symmetry.h"

using std::set;
using std::pair;
using std::swap;

const vec3d A3(1,1,1);               //3-fold axis for T, O and I symmetry types
const vec3d A5(0, 1, (sqrt(5)+1)/2); //5-fold axis for I symmetry type 

t_set &t_set::product(const t_set& s1, const t_set &s2)
{
   clear();
   set<mat3d>::iterator i1, i2;
   for(i1=s1.trans.begin(); i1!=s1.trans.end(); i1++)
      for(i2=s2.trans.begin(); i2!=s2.trans.end(); i2++)
         add((*i1) * (*i2));
   return *this;
}

t_set &t_set::product_with(const t_set& s1)
{
   set<mat3d> trans_orig = trans;
   clear();
   set<mat3d>::iterator i0, i1;
   for(i0=trans_orig.begin(); i0!=trans_orig.end(); i0++)
      for(i1=s1.trans.begin(); i1!=s1.trans.end(); i1++)
         add((*i0) * (*i1));
   return *this;
}

t_set &t_set::conjugate(const mat3d &m)
{
   set<mat3d> conj;
   mat3d inv_m = mat3d::inverse(m);
   for(set<mat3d>::iterator si=trans.begin(); si!=trans.end(); si++)
      conj.insert(m *(*si)*inv_m);
   trans = conj;
   return *this;
}

t_set &t_set::intersection(const t_set &s1, const t_set &s2)
{
   clear();
   set_intersection(s1.trans.begin(), s1.trans.end(),
                    s2.trans.begin(), s2.trans.end(),
                    inserter(trans, trans.end()) );
   return *this;
}

t_set &t_set::subtract(const t_set &ts)
{
   set<mat3d> diff;
   set_difference(trans.begin(), trans.end(),
                  ts.trans.begin(), ts.trans.end(),
                  inserter(diff, diff.begin()) );
   trans = diff;
   return *this;
}


t_set t_set::lcoset(const mat3d &m) const
{
   t_set coset;
   for(set<mat3d>::iterator si=trans.begin(); si!=trans.end(); si++)
      coset.add(m*(*si));
   return coset;
}

t_set operator *(const t_set &s1, const t_set &s2)
{ return t_set().product(s1, s2); }

t_set &operator *=(t_set &s1, const t_set &s2)
{ return s1.product_with(s2); }

t_set operator +(const t_set &s, const mat3d &m)
{ return t_set(s).add(m); }

t_set &operator +=(t_set &s, const mat3d &m)
{ return s.add(m); }
 
int compare(const t_set &t0, const t_set &t1)
{
   if(t0.size()>t1.size())
      return 1;
   else if(t0.size()<t1.size())
      return -1;

   // same size so test matrices one at a time
   set<mat3d>::const_iterator i0 = t0.begin();
   set<mat3d>::const_iterator i1 = t1.begin();
   for( ; i0!=t0.end(); ++i0, ++i1) {
      int cmp = compare(*i0, *i1);
      if(cmp)
         return cmp;
   }

   return 0;
}




t_set &t_set::min_set(const t_set &tr_whole, const t_set &tr_part,
      const mat3d &pos)
{
   t_set whole = tr_whole;
   t_set part = tr_part;
   part.conjugate(pos);
   t_set inter;
   inter.intersection(part, whole);
   inter += mat3d();                      // must include unit!
   while(whole.is_set()) {
      mat3d tr = *whole.trans.begin();    // select a transf from init list
      add(tr);                            // and add it to the final list
      t_set coset = inter.lcoset(tr);     // find equivalent transformations
      whole.subtract(coset);              // remove them
   }
   return *this;
}


///Class for working with Schoenflies notation
class sch_gen: public t_set
{
   private:
      sch_gen &operator *(const t_set &s) { product_with(s); return *this; }
      sch_gen &operator +(const mat3d &m) { add(m);  return *this; }
      
      ///Set up mirror transformation group.
      /**\param norm The normal for the mirror.
       *\return reference to this object with the transformations set. */
      sch_gen &refl(const vec3d &norm) { return unit() + mat3d::refl(norm); }
    
   public:
     
      ///Constructor
      /**\param t initialise with these transformations. */
      sch_gen(t_set t=t_set()): t_set(t) {};

      ///Set up unit symmetry transformation group.
      /**No relevant alignment.
       *\return reference to this object with the transformation set. */
      sch_gen &unit(){ add(mat3d()); return *this; }

      ///Set up horizontal mirror transformation group.
      /**Mirror normal (0,0,1).
       *\return reference to this object with the transformations set. */
      sch_gen &h_refl() { return refl(vec3d::Z); }
      
      ///Set up vertical mirror transformation group.
      /**Mirror normal (0,1,0).
       *\return reference to this object with the transformations set. */
      sch_gen &v_refl() { return refl(vec3d::Y); }

      ///Set up vertical mirror transformation group.
      /**Mirror normal (0,1,0) rotated PI/2n radians around (0,0,1).
       *\param n used for angle to rotate mirror, as PI/n radians.
       *\return reference to this object with the transformations set. */
      sch_gen &v_refl(int n)
         { return refl(mat3d::rot(vec3d::Z,-0.5*M_PI/n)*vec3d::Y); }

      ///Set up dihedral symmetry transformation group.
      /**Rotation axis (1,0,0).
       *\return reference to this object with the transformations set. */
      sch_gen &C2()
         { C(2); conjugate(mat3d::rot(vec3d::Z, vec3d::X)); return *this; }

      ///Set up Cs transformation group.
      /**Mirror normal in direction (0,0,1).
       *\return reference to this object with the transformations set. */
      sch_gen &Cs() { return h_refl(); }

      ///Set up Ci transformation group.
      /**No relevant alignment.
       *\return reference to this object with the transformations set. */
      sch_gen &Ci() { return unit() + mat3d::inversion(); }

      ///Set up C symmetry transformation group.
      /**Principal axis (0,0,1).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &C(int n);

      ///Set up Cv symmetry transformation group.
      /**Principal axis (0,0,1). Mirror normal (0,1,0).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &Cv(int n) { return C(n) * sch_gen().v_refl(); }

      ///Set up Ch symmetry transformation group.
      /**Principal axis (0,0,1). Mirror normal in direction (0,0,1).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &Ch(int n) { return C(n) * sch_gen().h_refl(); }

      ///Set up D symmetry transformation group.
      /**Principal axis (0,0,1). Dihedral axis in direction (1,0,0).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &D(int n) { return C(n) * sch_gen().C2(); }

      ///Set up Dv symmetry transformation group.
      /**Principal axis (0,0,1). Dihedral axis (1,0,0).
       * Vertical mirror normal (0,1,0) rotated PI/2n radians around (0,0,1).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &Dv(int n) { return D(n) * sch_gen().v_refl(n); }

      ///Set up Dh symmetry transformation group.
      /**Principal axis (0,0,1). Dihedral axis (1,0,0).
       * Horizontal mirror normal (0,0,1).
       *\param n principal axis is n-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &Dh(int n) { return D(n) * sch_gen().h_refl(); }

      ///Set up Sn symmetry transformation group.
      /**Principal axis (0,0,1). Horizontal mirror normal (0,0,1).
       *\param n principal axis is rotational n/2-fold.
       *\return reference to this object with the transformations set. */
      sch_gen &S(int n)
         { return C(n/2) * (sch_gen().unit() +
               mat3d::refl(vec3d::Z)*mat3d::rot(vec3d::Z,2*M_PI/n)); }

      ///Set up T symmetry transformation group.
      /**3-fold axes (1,1,1), (1,-1,-1).
       *\return reference to this object with the transformations set. */
      sch_gen &T()
         { return D(2) * (sch_gen().C(3).conjugate(mat3d::rot(vec3d::Z, A3))); }

      ///Set up Td symmetry transformation group.
      /**3-fold axes (1,1,1), (1,-1,-1).
       *\return reference to this object with the transformations set. */
      sch_gen &Td()
         { return T() * sch_gen().refl(vec3d::Y+vec3d::X); }

      ///Set up Th symmetry transformation group.
      /**3-fold axes (1,1,1), (1,-1,-1).
       *\return reference to this object with the transformations set. */
      sch_gen &Th() { return T() * sch_gen().Ci(); }

      ///Set up O symmetry transformation group.
      /**4-fold axes (1,0,0), (0,1,0).
       *\return reference to this object with the transformations set. */
      sch_gen &O()
         { return T() * (sch_gen().unit() + mat3d::rot(vec3d::X, M_PI/2)); }

      ///Set up Oh symmetry transformation group.
      /**4-fold axes (1,0,0), (0,1,0).
       *\return reference to this object with the transformations set. */
      sch_gen &Oh() { return O() * sch_gen().Ci(); return *this; }

      ///Set up I symmetry transformation group.
      /**5-fold axes (0,1,phi), (0,1,-phi).
       *\return reference to this object with the transformations set. */
      sch_gen &I()
         { return T() * sch_gen().C(5).conjugate(mat3d::rot(vec3d::Z, A5)); }

      ///Set up Ih symmetry transformation group.
      /**5-fold axes (0,1,phi), (0,1,-phi).
       *\return reference to this object with the transformations set. */
      sch_gen &Ih() { return I() * sch_gen().Ci(); }
};

sch_gen &sch_gen::C(int n)
{
   clear();
   unit();
   for(int i=1; i<n; i++)
      add(mat3d::rot(vec3d::Z, 2*M_PI*i/n));
   return *this;
}



void iso_type::normalise()
{
   if(axis.is_set()) {
      if( (axis[0]<-sym_eps) ||
          (axis[0]<sym_eps && axis[1]<-sym_eps) ||
          (axis[0]<sym_eps && axis[1]<sym_eps && axis[2]<-sym_eps) ) {
         axis = -axis;
         ang = -ang;
      }
      if(fabs(axis[0])>sym_eps || fabs(axis[1])>sym_eps ||fabs(axis[2])>sym_eps)
         axis.to_unit();
      else
         axis = vec3d();
   }

   ang = fmod(ang, 2*M_PI);
   if(ang<0) {
      if(ang>-sym_eps)
         ang = 0;
      else
         ang += 2*M_PI;
   }
}


bool iso_type::is_isometry(mat3d m) const
{
   vec3d rows[3];
   rows[0] = vec3d(m[0], m[1], m[2]);
   rows[1] = vec3d(m[4], m[5], m[6]);
   rows[2] = vec3d(m[8], m[9], m[10]);
   for(int i=0; i<3; i++)
      for(int j=i; j<3; j++) {
         double dot = vdot(rows[i], rows[j]);
         /// dot should equal 1 if i==j, otherwise 0
         if(dot<(i==j)-sym_eps || dot>(i==j)+sym_eps)
            return false;
      }

   return true;
}

iso_type &iso_type::init(mat3d m)
{
   if(!is_isometry(m)) {
      rot_type = rt_none;
      return *this;
   }

   transl = m.get_transl();
   m[3] = m[7] = m[11] = 0;  // zero translation column
   double det = m.det();
   if(det<0)
      m *= mat3d::inversion();

   vec4d quat = m.get_quaternion();
   double cos_a = quat[3];
   ang = 2*acos(safe_for_trig(cos_a));

   if(sqrt(1-cos_a*cos_a)>epsilon)
      axis = vec3d(quat[0], quat[1], quat[2]);
   else
      axis.unset();

   if(det>0) {
      if(fabs(ang)<sym_eps)
         rot_type = rt_unit;
      else
         rot_type = rt_rot;
   }
   else {
      if(fabs(ang)<sym_eps)
         rot_type = rt_inv;
      else if(fabs(fabs(ang)-M_PI)<sym_eps)
         rot_type = rt_refl;
      else
         rot_type = rt_rot_refl;
      ang = M_PI+ang;
      }

      normalise();
      return *this;
}
     
void iso_type::dump() const
{
   const char *et[] = { "none", "unit", "rot", "inv", "refl", "rot_refl"};
   fprintf(stderr, "rot_type=%s", et[rot_type]);
   if(rot_type==rt_rot || rot_type==rt_refl || rot_type==rt_rot_refl) {
      fprintf(stderr, ", axis=(%-.3f, %-.3f, %-.3f)",
            axis[0], axis[1], axis[2]);
      if(rot_type==rt_rot || rot_type==rt_rot_refl)
         fprintf(stderr, ", ang=%3.0f", ang*180/M_PI);
   }
   if(rot_type != rt_none)
      fprintf(stderr, ", transl=(%-.3f, %-.3f, %-.3f)",
            transl[0], transl[1], transl[2]);
   fprintf(stderr, "\n");
}



sch_axis::sch_axis(const mat3d &m)
{
   iso_type rot(m);
   double ang = rot.get_ang();
   if(fabs(ang)<sym_eps)
      nfold = 1;
   else {
      long tmp;
      double2rational(ang/(2*M_PI), tmp, nfold, sym_eps);
   }
   axis = rot.get_axis();
   //map index {rt_none=0, rt_unit, rt_rot, rt_inv, rt_refl, rt_rot_refl}
   int sym_map[] = { sch_sym::C1, sch_sym::C1, sch_sym::C,
                     sch_sym::Ci, sch_sym::Cs, sch_sym::S };
   sym_type = sym_map[rot.get_rot_type()];
}
    

void sch_axis::dump() const
{
   fprintf(stderr, "axis=%s", sch_sym(sym_type, nfold).get_symbol().c_str());
   if(axis.is_set())
      fprintf(stderr, ": a=(% .3f,% .3f,% .3f)", axis[0], axis[1], axis[2]);
   if(perp.is_set())
      fprintf(stderr, ", p=(% .3f,% .3f,% .3f)", perp[0], perp[1], perp[2]);
   fprintf(stderr, "\n");
}


bool sch_axis::operator <(const sch_axis &s) const
{
   int cmp = compare(axis, s.axis, sym_eps);
   if(cmp<0)
      return true;
   else if(cmp==0) {
      if(sym_type < s.sym_type)
         return true;
      else if(sym_type==s.sym_type) {
         if(nfold < s.nfold)
            return true;
         else if(nfold==s.nfold) {
            int cmp = compare(perp, s.perp, sym_eps);
            if(cmp<0)
               return true;
         }
      }
   }
   return false;
}




//find symmetry from rotational axes
void sch_sym::find_full_sym_type(const set<sch_axis> &full_sym)
{
   // find two greatest n-fold axes
   sch_axis max_fold1, max_fold2;
   bool has_dv=false, has_dh=false;
   for(set<sch_axis>::iterator si=full_sym.begin(); si!=full_sym.end(); si++) {
      int nfold = si->get_nfold();
      if(nfold>max_fold1.get_nfold())
         max_fold1 = *si;
      else if (nfold>max_fold2.get_nfold())
         max_fold2 = *si;
      if(si->get_sym_type()==sch_sym::Dv)
         has_dv = true;
      if(si->get_sym_type()==sch_sym::Dh)
         has_dh = true;
   }
   nfold = max_fold1.get_nfold();

   // dihedral axes only
   if(max_fold1.get_nfold()==2) {
      sch_axis ax = max_fold1;
      for(set<sch_axis>::iterator si=full_sym.begin(); si!=full_sym.end(); si++)
         if(si->get_sym_type()==sch_sym::Dv) {
            ax = *si;
            break;
         }
      sym_type = ax.get_sym_type();
      to_std = mat3d::alignment(ax.get_axis(), ax.get_perp(),
               vec3d::Z, vec3d::X);
   }
   
   // principal axis
   else if(max_fold1.get_nfold()>2 && max_fold2.get_nfold()<=2) {
      sym_type = max_fold1.get_sym_type();
      vec3d perp = max_fold1.get_perp();
      if(perp.is_set())
         to_std = mat3d::alignment(max_fold1.get_axis(), max_fold1.get_perp(),
               vec3d::Z, vec3d::X);
      else
         to_std = mat3d::rot(max_fold1.get_axis(), vec3d::Z);
   }
   
   // tetrahedral
   else if(max_fold1.get_nfold()==3 || max_fold1.get_nfold()==6 ) {
      if(has_dh)
         sym_type = sch_sym::Th;
      else if(has_dv)
         sym_type = sch_sym::Td;
      else
         sym_type = sch_sym::T;
     
      vec3d A3b = mat3d::rot(vec3d::Z, M_PI)*A3;
      if(vdot(max_fold1.get_axis(), max_fold2.get_axis()) > 0)
         A3b *= -1;


      to_std = mat3d::alignment(max_fold1.get_axis(), max_fold2.get_axis(),
            A3, A3b);
   }
   
   // octahedral
   else if(max_fold1.get_nfold()==4) {
      if(has_dh)
         sym_type = sch_sym::Oh;
      else
         sym_type = sch_sym::O;
      
      to_std = mat3d::alignment(max_fold1.get_axis(), max_fold2.get_axis(),
                           vec3d::Z, vec3d::X);
   }
   
   // icosahedral
   else if(max_fold1.get_nfold()==5) {
      if(has_dh)
         sym_type = sch_sym::Ih;
      else
         sym_type = sch_sym::I;
      
      vec3d A5b = mat3d::rot(A3, 2*M_PI/3) * A5;
      if(vdot(max_fold1.get_axis(), max_fold2.get_axis()) < 0)
         A5b *= -1;
      
      to_std = mat3d::alignment(max_fold1.get_axis(), max_fold2.get_axis(),
                           A5, A5b);
   }

   // unknown symmetry
   else {
      //sym_type = sch_sym::unknown;
      //Don't allow to fail
      sym_type = sch_sym::C1;
      to_std = mat3d();
   }
}

static bool operator <(const vec3d &v1, const vec3d &v2)
{
   return compare(v1, v2, sym_eps)<0;
}

sch_sym::sch_sym(const t_set &ts): sym_type(unknown), nfold(1), to_std(mat3d())
{
   if(!ts.size())
      return;

   // Translate the set of symmetry elements to fix the origin
   vec3d fixed_pt(0,0,0);
   for(t_set::const_iterator ti=ts.begin(); ti!=ts.end(); ++ti)
      fixed_pt += vec3d((*ti)[3], (*ti)[7], (*ti)[11]);
   fixed_pt /= ts.size();  // centroid of points where origin is sent
   
   mat3d transl = mat3d::transl(-fixed_pt);
   t_set o_ts = ts;
   o_ts.conjugate(transl);

  
   for(t_set::const_iterator ti=o_ts.begin(); ti!=o_ts.end(); ++ti)
      axes.insert(sch_axis(*ti));

   bool inv = false;
   vector<vec3d> refls, dihs;
   for(set<sch_axis>::const_iterator si=axes.begin(); si!=axes.end(); si++) {
      if(si->get_sym_type()==sch_sym::Cs)
         refls.push_back(si->get_axis());
      else if(si->get_sym_type()==sch_sym::C && si->get_nfold()==2)
         dihs.push_back(si->get_axis());
      else if(si->get_sym_type()==sch_sym::Ci)
            inv = true;
   }
  
   map<vec3d, vector<sch_axis> > v2ax;
   for(set<sch_axis>::const_iterator si=axes.begin(); si!=axes.end(); si++) {
      const vec3d &axis = si->get_axis();
      if(axis.is_set())
         v2ax[axis].push_back(*si);
   }

   if(v2ax.size()==0) { // no axes
      axes.clear();
      if(inv)
         sym_type = sch_sym::Ci;
      else
         sym_type = sch_sym::C1;
   }
   else if(v2ax.size()==1 && v2ax.begin()->second.size()==1 &&
          v2ax.begin()->second[0].get_sym_type()==sch_sym::Cs) { //one refl
      axes.clear();
      to_std = mat3d::rot(v2ax.begin()->second[0].get_axis(), vec3d::Z);
      sym_type = sch_sym::Cs;
   }
   else {// remaining possibilities have rotational axis
      set<sch_axis> full_sym;
      map<vec3d, vector<sch_axis> >::iterator mi;
      for(mi=v2ax.begin(); mi!=v2ax.end(); mi++) {
         vec3d axis = mi->first;
         int nfold=0;
         int sym = -1;
         int h_refl = false;
         for(vector<sch_axis>::iterator vi = mi->second.begin();
               vi!=mi->second.end(); vi++) {
            if(vi->get_sym_type()==sch_sym::Cs)
               h_refl = true;
            if(vi->get_nfold()>nfold) {
               nfold = vi->get_nfold();
               sym = vi->get_sym_type();
            }
         }

         sch_axis sym_ax;
         sym_ax.set_axis(axis);
         if(!h_refl) {              // no horz reflection -> S or C
            sym_ax.set_sym_type(sym);
            sym_ax.set_nfold(nfold);
         }
         else {                     // horz refl -> Ch
            if(sym_ax.get_sym_type()==sch_sym::S)
               sym_ax.set_nfold(nfold/2);
            else
               sym_ax.set_nfold(nfold);
            sym_ax.set_sym_type(sch_sym::Ch);
         }

         // vertical reflection
         for(vector<vec3d>::iterator vi=refls.begin(); vi!=refls.end(); ++vi) {
            if(fabs(vdot(*vi, axis))<sym_eps) {
               if(sym_ax.get_sym_type()==sch_sym::C) {
                  sym_ax.set_sym_type(sch_sym::Cv);
                  sym_ax.set_perp(mat3d::rot(axis, M_PI/2)*(*vi));
               }
               else if(sym_ax.get_sym_type()==sch_sym::Ch) {
                  sym_ax.set_sym_type(sch_sym::Dh);
                  sym_ax.set_perp(mat3d::rot(axis, M_PI/2+M_PI/nfold)*(*vi));
               }
               else {  // sch_sym::S
                  sym_ax.set_sym_type(sch_sym::Dv);
                  sym_ax.set_nfold(nfold/2);
                  sym_ax.set_perp(mat3d::rot(axis,
                           M_PI/2-0.5*M_PI/sym_ax.get_nfold())*(*vi));
               }
               break;
            }
         }

         // dihedral axis
         if(sym_ax.get_sym_type()!=sch_sym::Dv && sym_ax.get_sym_type()!=sch_sym::Dh) {
            for(vector<vec3d>::iterator vi=dihs.begin(); vi!=dihs.end(); ++vi) {
               if(fabs(vdot(*vi, axis))<sym_eps) {
                  if(sym_ax.get_sym_type()==sch_sym::C) {
                     sym_ax.set_perp(*vi);
                     sym_ax.set_sym_type(sch_sym::D);
                  }
                  else { // sym==sch_sym::Cv
                     sym_ax.set_sym_type(sch_sym::Dv);
                     sym_ax.set_perp(mat3d::rot(axis,
                              -0.5*M_PI/sym_ax.get_nfold())*(sym_ax.get_perp()));
                  }
                  break;
               }
            }
         }
         if(nfold>1) {               // has rotation axis (not just refl)
            full_sym.insert(sym_ax);
            //sym_ax.dump();
         }
      }

      axes = full_sym;
      find_full_sym_type(full_sym);
   }
  
   to_std *= transl;  // first move the fixed point to the origin;
}



static inline int edge_seen(map<pair<int, int>, unsigned char> &e_seen, int v0, int v1, bool update)
{
   unsigned char dir = 1;
   if(v0>v1) {
      swap(v0, v1);
      dir = 2;
   }
   pair<int, int> edge(v0, v1);
   map<pair<int, int>, unsigned char>::iterator ei = e_seen.find(edge);
   if(ei == e_seen.end()) {  // not traversed in either direction
      if(update)
         e_seen[edge] = dir; // mark as traversed in this direction
      return 0;              // 0: unseen
   }
   if(ei->second & dir)
      return 2;              // 2: seen and already traversed in this direction
   if(update)
      ei->second |= dir;     // mark as traversed in this direction
   return 1;                 // 1: seen but not traversed in this direction
}

static inline int edge_check(map<pair<int, int>, unsigned char> &e_seen, int v0,int v1)
{ return edge_seen(e_seen, v0, v1, false); }

static inline int edge_mark(map<pair<int, int>, unsigned char> &e_seen, int v0, int v1)
{ return edge_seen(e_seen, v0, v1, true); }
   

// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.30.6536
// Symmetries of Polyhedra: Detection and Applications
// by X. Y. Jiang, H. Bunke
// ftp://ftp.iam.unibe.ch/pub/TechReports/1994/iam-94-012.ps.gz 
// 3.1.1 The algorithm of Jiang & Bunke

static int find_path(vector<int> &path, vector<int> &v_code,
      const vector<int> &edge, const vector<vector<int> > &v_cons,
      const vector<int> *test_path=0, const vector<int> *test_v_code=0)
{
   path.clear();
   v_code = vector<int>(v_cons.size(), -1);
   map<pair<int, int>, unsigned char> e_seen;
   int v_cnt = 0;
   int v_cur = edge[0];
   path.push_back(v_cur);
   v_code[v_cur] = v_cnt++;
   int v_next = edge[1];
   while(true) {
      int seen = edge_mark(e_seen, v_cur, v_next);
      int v_new = -1;
      vector<int>::const_iterator vi, vi_orig;
      if(v_code[v_next] < 0) { // new vertex, exit from "right"
         vi = find(v_cons[v_next].begin(), v_cons[v_next].end(), v_cur);
         if(++vi ==  v_cons[v_next].end())
            vi = v_cons[v_next].begin();
         v_code[v_next] = v_cnt++;
         v_new = *vi;
      }
      else { // previously seen vertex
         if(seen) { // previously seen edge, exit from "right"
            vi_orig = find(v_cons[v_next].begin(), v_cons[v_next].end(), v_cur);
            for(vi=vi_orig+1; vi!=vi_orig; ++vi) {  // check for an edge out
               if(vi ==  v_cons[v_next].end()) {
                  vi = v_cons[v_next].begin();      // wrap to beginning
                  if(vi == vi_orig)
                     break;
               }
               if(edge_check(e_seen, v_next, *vi)<2) {
                  v_new = *vi;
                  break;
               }
            }
            if(vi == vi_orig) // finished, didn't find an edge to leave from
               return 1;
         }
         else { // new edge, go back
            v_new = v_cur;
         }
      }

      // stop if the coded path is different to the test coded path
      if(test_path && v_code[v_next]!=(*test_v_code)[(*test_path)[path.size()]])
         return 0;
      
      path.push_back(v_next);
      v_cur = v_next;
      v_next = v_new;
   }
}  


static void update_equiv_elems(vector<map<int, set<int> > > &equiv_elems,
      const vector<map<int, set<int> > > &new_equivs, int cnts[3])
{
   for(int i=0; i<3; i++) {
      map<int, set<int> >::const_iterator mi_new;
      for(mi_new=new_equivs[i].begin(); mi_new!=new_equivs[i].end(); ++mi_new) {
         int to = *mi_new->second.begin();
         if(to>=cnts[i])
            to -= cnts[i];
         for(set<int>::iterator si=mi_new->second.begin();
               si!=mi_new->second.end(); ++si) {
            int from = (*si<cnts[i]) ? *si : *si-cnts[i];
            equiv_elems[i][to].insert(from);
         }
      }
   }
}



static void equiv_elems_to_sets(vector<vector<set<int> > > &equiv_sets,
      vector<map<int, set<int> > > &equiv_elems,
      vector<map<int, set<int> > > &orig_equivs)
{
   equiv_sets.clear();
   equiv_sets.resize(3);
   for(int i=0; i<3; i++) {
      // map each index to the to minimum index it is equivalent to
      map<int, int> idx_to_min_idx;
      map<int, set<int> >::const_iterator mi;
      for(mi=equiv_elems[i].begin(); mi!=equiv_elems[i].end(); ++mi) {
         set<int>::const_iterator si = mi->second.begin();
         int min_idx = *si;
         for( ;si!=mi->second.end(); ++si) {
            map<int, int>::iterator idx_i = idx_to_min_idx.find(*si);
            if(idx_i == idx_to_min_idx.end())
               idx_to_min_idx[*si] = min_idx;
            else if(idx_i->second > min_idx)
               idx_i->second = min_idx;
         }
      }

      // map each minimum index to its set of equivalent indexes
      map<int, set<int> > min_idx_to_idxs;
      map<int, int>::const_iterator idx_i;
      for(idx_i=idx_to_min_idx.begin(); idx_i!=idx_to_min_idx.end(); ++idx_i)
         min_idx_to_idxs[idx_i->second].insert(idx_i->first);

      // Copy sets to the return argument
      map<int, set<int> >::const_iterator midx_i;
      for(midx_i=min_idx_to_idxs.begin(); midx_i!=min_idx_to_idxs.end();
            ++midx_i) {
         equiv_sets[i].push_back(set<int>());
         for(set<int>::const_iterator si=midx_i->second.begin();
               si!=midx_i->second.end(); ++si)
               for(set<int>::const_iterator si2=equiv_elems[i][*si].begin();
                     si2!=equiv_elems[i][*si].end(); ++si2)
                  equiv_sets[i].back().insert(
                      orig_equivs[i][*si2].begin(), orig_equivs[i][*si2].end());
      }
   }

}




static bool is_sym(const geom_if &test_geom, const geom_if &geom,
      const vector<int> &test_v_code, const vector<int> &v_code,
      bool orient, mat3d &trans, vector<map<int, set<int> > > &new_equivs)
{
   int v_sz = test_geom.verts().size();
   // code to vertex idx for this sym
   vector<int> c2v_map(v_sz);
   for(int v=0; v<v_sz; v++)
      c2v_map[v_code[v]] = v;
   vector<int> v_map(v_sz);
   for(int v=0; v<v_sz; v++)
      v_map[v] = c2v_map[test_v_code[v]];

   vector<vec3d> t_pts(3), pts(3);
   for(int i=0; i<2; i++) {
      t_pts[i] = test_geom.verts(i);
      pts[i]   = test_geom.verts(v_map[i]);
   }
   // Choose a third point that is not colinear with the first two
   for(int i=2; i<v_sz; i++) {
      t_pts[2] = test_geom.verts(i);
      pts[2]   = test_geom.verts(v_map[i]);
      vec3d norm = vcross(pts[1]-pts[0], pts[2]-pts[0]);
      if(norm.mag2()>(epsilon*epsilon))
         break;
   }

   if(orient)
      transform(pts, mat3d::inversion());
   trans = mat3d::alignment(t_pts, pts);
   if(orient)
      trans = mat3d::inversion()*trans;
   col_geom_v s_geom = geom;
   s_geom.transform(trans);

   bool is_congruent = check_congruence(geom, s_geom, &new_equivs, sym_eps);
   return is_congruent;
}

static void set_equiv_elems_identity(const geom_if &geom,
      vector<vector<set<int> > > *equiv_sets)
{
   int cnts[3] = { geom.verts().size(),
                   geom.edges().size(),
                   geom.faces().size() };
   equiv_sets->clear();
   equiv_sets->resize(3);
   for(int i=0; i<3; i++) {
      (*equiv_sets)[i].resize(cnts[i]);
      for(int j=0; j<cnts[i]; j++)
         (*equiv_sets)[i][j].insert(j);
   }
}

static int find_syms(const geom_if &geom, t_set &ts,
      vector<vector<set<int> > > *equiv_sets)
{
   ts.clear();
   
   col_geom_v merged_geom = geom;
   vector<map<int, set<int> > > orig_equivs;
   sort_merge_elems(merged_geom, "vef", &orig_equivs, epsilon);
   col_geom_v test_geom = merged_geom;

   geom_info inf(merged_geom);
   //if(inf.num_parts()>1 || !inf.is_orientable()) // for octahemioctaheron=Td
   const int dim = test_geom.set_hull(msg_str("-A%.15f", 1.0-sym_eps));
   if(dim<2)  // contains an infinite axis, can't currently handle this
      return 0;

   geom_info g_inf(test_geom);
   const vector<vector<int> > &edges = g_inf.get_impl_edges();
   const vector<vector<int> > &v_cons = g_inf.get_vert_cons();
   
   vector<vector<int> > r_cons = v_cons;
   for(vector<vector<int> >::iterator vi=r_cons.begin(); vi!=r_cons.end(); vi++)
      reverse(vi->begin(), vi->end());
   const vector<vector<int> > *cons[] = {&v_cons, &r_cons};
   
   vector<map<int, set<int> > > equiv_elems(3);
   int cnts[3] = { merged_geom.verts().size(),
                   merged_geom.edges().size(),
                   merged_geom.faces().size() };
   vector<int> test_path, path;
   vector<int> test_v_code, v_code;
   find_path(test_path, test_v_code, *edges.begin(), v_cons);
   vector<vector<int> >::const_iterator ei;
   for(ei=edges.begin(); ei!=edges.end(); ei++) {
      vector<int> edge = *ei;
      for(int i=0; i<2; i++) {
         if(i)
            swap(edge[0], edge[1]);
         for(int orient=0; orient<2; orient++) {
            if(find_path(path, v_code, edge, *cons[orient],
                     &test_path, &test_v_code)) {
               mat3d trans;
               vector<map<int, set<int> > > new_equivs;
               if(is_sym(test_geom, merged_geom, test_v_code,v_code, orient,
                        trans, new_equivs)) {
                  ts.add(trans);
                  if(equiv_sets)
                     update_equiv_elems(equiv_elems, new_equivs, cnts);
               }
            }
         }
      }
   }

   if(equiv_sets)
      equiv_elems_to_sets(*equiv_sets, equiv_elems, orig_equivs);

   // Don't allow to fail
   if(ts.size()==0) {
      ts.add(mat3d());
      if(equiv_sets)
         set_equiv_elems_identity(geom, equiv_sets);
   }

   return 1;
}



const char *type_str[sch_sym::Ih+1] = {
      "Unknown", "C1", "Ci", "Cs", "C", "Cv", "Ch", "D", "Dv", "Dh", "S",
      "T", "Td", "Th", "O", "Oh", "I", "Ih" };


sch_sym::sch_sym(const geom_if &geom,
      vector<vector<set<int> > > *equiv_sets)
{
   init(geom, equiv_sets);
}

sch_sym::sch_sym(int type, int n, const mat3d &pos, char *errmsg)
{
   init(type, n, pos, errmsg);
}

sch_sym::sch_sym(const string &name, const mat3d &pos, char *errmsg)
{
   init(name, pos, errmsg);
}
      
sch_sym::sch_sym(const sch_axis &sym_axis, const vec3d &cent)
{
   init(sym_axis, cent);
}


bool sch_sym::init(const geom_if &geom,
      vector<vector<set<int> > > *equiv_sets)
{
   sym_type = unknown;
   t_set ts;
   find_syms(geom, ts, equiv_sets);
   *this = sch_sym(ts);
   return sym_type!=unknown;
}


bool sch_sym::init(int type, int n, const mat3d &pos, char *errmsg)
{
   sym_type = unknown;
   nfold = 0;
   axes.clear();
   mirrors.clear();
   sub_syms.clear();
   autos = sch_sym_autos();

   if((type<C1 || type>Ih)) {
      if(errmsg)
         strcpy(errmsg, "unknown symmetry type");
      return false;
   }
   else if((type<C || type>S) && n) {
      if(errmsg)
         strcpy(errmsg, "symmetry type doesn't take a number");
      return false;
   }
   else if((type>=C && type<=S) && !n) {
      if(errmsg)
         strcpy(errmsg, "symmetry type needs a number");
      return false;
   }
   else if (type==S && !is_even(n)) {
      if(errmsg)
         strcpy(errmsg, "symmetry type S needs an even number");
      return false;
   }

   to_std = pos;
   sym_type = type;
   if(sym_type>=C || sym_type<=S)   // principal axis
      nfold = n;

   //Allow for alternative descriptions of "lesser" symmetries
   if( (nfold<2 && sym_type>=C && sym_type<=Dh) ||
          (nfold==2 && sym_type==S) ) {
      sch_sym tmp(get_trans());
      *this = tmp;
   }


   // Maybe normalise symmetry, have to adjust to_std with this
   
   return true;

}

   
bool sch_sym::init(const string &name, const mat3d &pos, char *errmsg)
{
   sym_type = unknown;
   nfold = 0;

   if(name.length() > MSG_SZ) {
      if(errmsg)
         strcpy(errmsg, "too many characters in symmetry type name");
      return false;
   }

   char c_num[MSG_SZ] = "";
   char c_name[MSG_SZ] = "";
   char numbers[] = "0123456789";
   size_t pos1 = strcspn(name.c_str(), numbers);
   size_t pos2 = pos1;
   if(pos1 < name.length())
      pos2 += strspn(name.c_str()+pos1, numbers);

   strncpy(c_name, name.c_str(), pos1);
   strcat(c_name, name.c_str()+pos2);
   strncpy(c_num, name.c_str()+pos1, pos2-pos1);

   int fold = atoi(c_num);
   map<string, int> sym_names;
   for(int i=1; i<Ih+1; i++)
      sym_names[type_str[i]] = i;
   
   map<string, int>::iterator sni = sym_names.find(c_name);
   int type = (sni != sym_names.end()) ? sni->second : -1;

   return init(type, fold, pos, errmsg);
}

void sch_sym::init(const sch_axis &sym_axis, const vec3d &cent)
{
   axes.clear();
   mirrors.clear();
   sub_syms.clear();
   autos = sch_sym_autos();

   sym_type = sym_axis.get_sym_type();
   nfold = sym_axis.get_nfold();
   to_std = mat3d::transl(cent) *
            mat3d::alignment(vec3d::X, vec3d::Z,
                             sym_axis.get_axis(), sym_axis.get_perp());
   
}


string sch_sym::get_symbol() const
{
   if(sym_type<C || sym_type>S)
      return type_str[sym_type];
   else
      return string() + type_str[sym_type][0] + itostr(nfold) +
         (type_str[sym_type]+1);
}


t_set &sch_sym::get_trans(t_set &ts) const
{
   ts.clear();
   if(sym_type>=C1 && sym_type<=Ih) {
      if(sym_type<C || sym_type>S) {
         typedef sch_gen &(sch_gen::* PF_SYM)();
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

         map<int, PF_SYM>::iterator soi = sym_names_other.find(sym_type);
         ts = (sch_gen().*soi->second)();
      }
      else {
         typedef sch_gen &(sch_gen::* PF_SYMN)(int n);
         map<int, PF_SYMN> sym_names_n;
         sym_names_n[C] = &sch_gen::C;
         sym_names_n[Cv] = &sch_gen::Cv;
         sym_names_n[Ch] = &sch_gen::Ch;
         sym_names_n[D] = &sch_gen::D;
         sym_names_n[Dv] = &sch_gen::Dv;
         sym_names_n[Dh] = &sch_gen::Dh;
         sym_names_n[S] = &sch_gen::S;

         map<int, PF_SYMN>::iterator sni = sym_names_n.find(sym_type);
         ts = (sch_gen().*sni->second)(nfold);
      }
   }

   ts.conjugate(mat3d::inverse(to_std));
   return ts;
}


t_set sch_sym::get_trans() const
{ 
   t_set ts; return get_trans(ts);
}
   

const set<sch_axis> &sch_sym::get_axes() const
{
   if(axes.size()==0 && sym_type!=unknown) {
      t_set ts;
      get_trans(ts);
      axes = sch_sym(ts).axes;
   }
   return axes;
}

const set<vec3d> &sch_sym::get_mirrors() const
{
   if(mirrors.size()==0 &&
         ( sym_type==Cs || sym_type==Cv || sym_type==Ch ||
           sym_type==Dv || sym_type==Dh || sym_type==Td ||
           sym_type==Th || sym_type==Oh || sym_type==Ih ) ) {
      t_set ts;
      get_trans(ts);
      for(t_set::const_iterator ti=ts.begin(); ti!=ts.end(); ++ti) {
         sch_axis ax = *ti;
         if(ax.get_sym_type() == Cs)
            mirrors.insert(ax.get_axis());
      }
   }
   return mirrors;
}

bool sch_sym::has_inversion_symmetry() const
{
   if(sym_type==Ih||sym_type==Oh ||sym_type==Th || sym_type==Ci ||
         (sym_type==Dv&&nfold%2) ||
         (sym_type==Dh&&nfold%2==0) ||
         (sym_type==S&&(nfold/2)%2) )
      return true;
   else
      return false;
}


void get_equiv_elems(const geom_if &geom, const t_set &ts,
      vector<vector<set<int> > > *equiv_sets)
{

   equiv_sets->clear();
   if(ts.size()<=1) {
      set_equiv_elems_identity(geom, equiv_sets);
      return;
   }

   col_geom_v merged_geom = geom;
   vector<map<int, set<int> > > orig_equivs;
   sort_merge_elems(merged_geom, "vef", &orig_equivs, epsilon);
   col_geom_v test_geom = merged_geom;

   vector<map<int, set<int> > > equiv_elems(3);
   int cnts[3] = { merged_geom.verts().size(),
                   merged_geom.edges().size(),
                   merged_geom.faces().size() };

   for(set<mat3d>::iterator si=ts.begin(); si!=ts.end(); si++) {
      col_geom_v trans_geom = merged_geom;
      trans_geom.transform(*si);
      vector<map<int, set<int> > > new_equivs;
      check_congruence(merged_geom, trans_geom, &new_equivs, sym_eps);
      update_equiv_elems(equiv_elems, new_equivs, cnts);
   }

   if(equiv_sets)
      equiv_elems_to_sets(*equiv_sets, equiv_elems, orig_equivs);
}


void sch_sym::add_sub_axes(const sch_sym &sub) const
{
   vector<int> factors;
   int fold=sub.get_nfold();
   for(int i=1; i<=fold/2; i++)
      if(fold%i==0)
         factors.push_back(fold/i);
   for(unsigned int i=0; i<factors.size(); i++) {
      int nfold = factors[i];
      sch_sym sub_ax(sub);
      sub_ax.set_nfold(nfold);
      switch(sub.get_sym_type()) {
         case C:
            sub_syms.insert(sub_ax);
            break;
         case Cv:
            sub_syms.insert(sub_ax);
            sub_ax.set_sym_type(C);
            sub_syms.insert(sub_ax);
            if(fold%2==0) {
               sub_ax.set_sym_type(Cv);
               sub_ax.set_to_std(mat3d::rot(0,0,M_PI/(fold))*sub.get_to_std());
               sub_syms.insert(sub_ax);
            }
            break;
         case Ch:
            sub_syms.insert(sub_ax);
            sub_ax.set_sym_type(C);
            sub_syms.insert(sub_ax);
            if(factors[i]>2 && factors[i]%2==0) {
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
            sub_ax.set_to_std(mat3d::rot(0,0,M_PI/(2*nfold))*sub.get_to_std());
            sub_syms.insert(sub_ax);
            if(nfold%2==0) {
               sub_ax.set_to_std(mat3d::rot(0,0,-M_PI/(2*nfold))*
                     sub.get_to_std());
               sub_syms.insert(sub_ax);
            }
            sub_syms.insert(sub_ax);
            sub_ax.set_nfold(2*nfold);
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
            if(nfold>2 && nfold%2==0) {
               sub_ax.set_sym_type(S);
               sub_syms.insert(sub_ax);
               sub_ax.set_nfold(nfold/2);
               sub_ax.set_sym_type(Dv);
               sub_syms.insert(sub_ax);
            }
            break;
         case S:
            if(nfold>2 && nfold%2==0) {
               sub_syms.insert(sub_ax);
               sub_ax.set_sym_type(C);
               sub_ax.set_nfold(nfold/2);
               sub_syms.insert(sub_ax);
            }
            break;

       }
   }
}


class same_sym_group {
   private:
      const sch_sym &sym;
   public:
      same_sym_group(sch_sym s): sym(s) {}
   
      bool operator () (const sch_sym &cmp) const {
         if(cmp.get_sym_type()==sym.get_sym_type()) {
            if(sym.get_sym_type()<sch_sym::C || sym.get_sym_type()>sch_sym::S)
               return true;
            else if(cmp.get_nfold()==sym.get_nfold())
               return true;
         }
         return false;
      }
};



sch_sym sch_sym::get_sub_sym(const string &sub, char *errmsg) const
{
   char errmsg2[MSG_SZ];
   sch_sym sub_sym, final_sym;
   
   char *sub_cpy = copy_str(sub.c_str());   
   bool valid = true;    // to allow free of sub_cpy in one place
   vector<char *> parts;
   split_line(sub_cpy, parts, ",");
   if(parts.size()>2) {
      strcpy(errmsg, "too many comma separated parts");
      valid = false;
   }
   
   if(valid) {
      if(parts.size()==0 || strncmp(parts[0], "full", strlen(parts[0]))==0) {
         sub_sym = *this;
      }
      else if(!sub_sym.init(parts[0], mat3d(), errmsg2)) {
         snprintf(errmsg, MSG_SZ, "sub-symmetry type: %s", errmsg2);
         valid = false;
      }
   }
    
   int sub_sym_conj = 0;
   if(valid && parts.size()>1) {
      if(!read_int(parts[1], &sub_sym_conj, errmsg2)) {
         snprintf(errmsg, MSG_SZ,"sub-symmetry conjugation number: %s",errmsg2);
         valid = false;
      }
   }

   if(valid) {
      final_sym = get_sub_sym(sub_sym, sub_sym_conj, errmsg2);
      if(final_sym.get_sym_type() == sch_sym::unknown) {
         snprintf(errmsg, MSG_SZ, "sub-symmetry: %s", errmsg2);
         valid = false;
      }
   }

   free(sub_cpy);
   return valid ? final_sym : sch_sym();
}


sch_sym sch_sym::get_sub_sym(const sch_sym &sub_sym, int conj_type,char *errmsg) const
{
   if(errmsg) {
      *errmsg = '\0';
   }
   if(conj_type<0) {
      if(errmsg)
         sprintf(errmsg, "conjugation type number cannot be negative");
      return sch_sym();
   }
   get_sub_syms();
   int conj = -1;
   set<sch_sym>::const_iterator si = sub_syms.begin();
   same_sym_group cmp(sub_sym);
   while((si = find_if(si, sub_syms.end(), cmp))!=sub_syms.end()) {
      ++conj;
      if(conj==conj_type)
         return *si;
      ++si;
   }
   if(conj<0 && errmsg)
      sprintf(errmsg, "%s is not a sub-symmetry of %s",
            sub_sym.get_symbol().c_str(), get_symbol().c_str());
   else if(si==sub_syms.end() && errmsg)
      sprintf(errmsg, "conjugation type too large for %s (last number: %d)",
            sub_sym.get_symbol().c_str(), conj);

   return sch_sym(sch_sym::unknown);
}

const set<sch_sym> &sch_sym::get_sub_syms() const
{
   sub_syms.clear();
   const vec3d axis = vec3d::Z;
   const vec3d perp = vec3d::X;
   int dih_type;
   if(sub_syms.size()==0 && sym_type!=unknown) {
      sch_sym sym = *this;
      sub_syms.insert(sym);
      if(sym_type!=C1) {    // Don't add C1 a second time
         sym.sym_type=C1;
         sub_syms.insert(sym);
      }
      if(has_inversion_symmetry()) {
         sym.sym_type=Ci;
         sub_syms.insert(sym);
      }
      switch(sym_type) {
         case Ih:
            sym.sym_type=I;
            sub_syms.insert(sym);
            sym.sym_type=Th;
            sub_syms.insert(sym);
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.sym_type=Ci;
            sub_syms.insert(sym);
            sym.sym_type=Cs;
            sub_syms.insert(sym);
            sym.init(Dv, 5, mat3d::alignment(A5, vec3d::Z,
                     axis, mat3d::rot(axis, M_PI/10)*vec3d::X) * to_std);
            add_sub_axes(sym);
            sym.init(Dv, 3, mat3d::alignment(A3, A5,
                     axis, mat3d::rot(axis, M_PI/6 )*vec3d::X) * to_std);
            add_sub_axes(sym);
            sym.init(Dh, 2, to_std);
            add_sub_axes(sym);
            break;

         case I:
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.init(D, 5, mat3d::alignment(A5, vec3d::Z,
                     axis, mat3d::rot(axis, M_PI/10)*vec3d::X) * to_std);
            add_sub_axes(sym);
            sym.init(D, 3, mat3d::alignment(A3, A5,
                     axis, mat3d::rot(axis, M_PI/6 )*vec3d::X) * to_std);
            add_sub_axes(sym);
            sym.init(D, 2, to_std);
            add_sub_axes(sym);
            break;
         
         case Oh:
            sym.sym_type=O;
            sub_syms.insert(sym);
            sym.sym_type=Td;
            sub_syms.insert(sym);
            sym.sym_type=Th;
            sub_syms.insert(sym);
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.sym_type=Ci;
            sub_syms.insert(sym);
            sym.sym_type=Cs;
            sub_syms.insert(sym);
            sym.init(Cs, 0, mat3d::rot(vec3d(0,1,1), axis) * to_std);
            add_sub_axes(sym);
            sym.init(Dh, 4, to_std);
            add_sub_axes(sym);
            sym.init(Dv, 3, mat3d::alignment(A3, vec3d(1,-1,0), axis, perp)
                  * to_std);
            add_sub_axes(sym);
            sym.init(Dh, 2, mat3d::rot(vec3d(0,1,1), axis) * to_std);
            add_sub_axes(sym);
            break;
            
         case O:
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.init(D, 4, to_std);
            add_sub_axes(sym);
            sym.init(D, 3, mat3d::alignment(A3, vec3d(1,-1,0), axis, perp)
                  * to_std);
            add_sub_axes(sym);
            sym.init(D, 2, mat3d::rot(vec3d(0,1,1), axis) * to_std);
            add_sub_axes(sym);
            break;

         case Th:
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.sym_type=Ci;
            sub_syms.insert(sym);
            sym.sym_type=Cs;
            sub_syms.insert(sym);
            sym.init(S, 6, mat3d::rot(A3, axis) * to_std);
            add_sub_axes(sym);
            sym.init(Dh, 2, to_std);
            add_sub_axes(sym);
            break;

         case Td:
            sym.sym_type=T;
            sub_syms.insert(sym);
            sym.init(Cs, 0, mat3d::rot(vec3d(0,1,1), axis) * to_std);
            sub_syms.insert(sym);
            sym.init(Cv, 3, mat3d::alignment(A3, vec3d(1,-1,-1), axis, perp)
                  * to_std);
            add_sub_axes(sym);
            sym.init(Dv, 2, to_std);
            add_sub_axes(sym);
            break;
         
         case T:
            sym.init(C, 3, mat3d::rot(A3, axis) * to_std);
            add_sub_axes(sym);
            sym.init(D, 2, to_std);
            add_sub_axes(sym);
            break;

         case Dh:
            sym.init(Cs, 0, to_std);     // horizontal mirror
            sub_syms.insert(sym);
            if(nfold%2==0) {       // nfold even: second axis and vert mirror
               dih_type = Dh;
               // vertical mirror through dihedral axis
               sym.init(Cs, 0, mat3d::rot(vec3d::X, -M_PI/2) *
                               mat3d::rot(axis, -M_PI/nfold) * to_std);
               sub_syms.insert(sym);
               sym.init(dih_type, 2, mat3d::rot(vec3d::Y, -M_PI/2) *
                                     mat3d::rot(axis, -M_PI/nfold) * to_std);
               add_sub_axes(sym);
               if(nfold>2) {
                  sym.init(Dh, nfold/2, mat3d::rot(axis, -M_PI/nfold) * to_std);
                  add_sub_axes(sym);
                  sym.init(Dv, nfold/2, mat3d::rot(axis, -M_PI/nfold) * to_std);
                  add_sub_axes(sym);
               }
            }
            else
               dih_type = Cv;

            // vertical mirror through dihedral axis
            sym.init(Cs, 0, mat3d::rot(vec3d::X, -M_PI/2) * to_std);
            sub_syms.insert(sym);
            sym.init(dih_type, 2, mat3d::rot(vec3d::Y, -M_PI/2) * to_std);
            add_sub_axes(sym);
            add_sub_axes(*this);
            break;

         case Dv:
            // vertical mirror between dihedral axes
            sym.init(Cs, 0, mat3d::rot(axis, -M_PI/(2*nfold)) *
                            mat3d::rot(vec3d::X, M_PI/2) * to_std);
            sub_syms.insert(sym);
            dih_type = nfold%2 ? Ch : D;
            sym.init(dih_type, 2, mat3d::rot(vec3d::Y, M_PI/2) * to_std);
            add_sub_axes(sym);
            sym.init(dih_type, 2, mat3d::rot(axis, M_PI/nfold) *
                            mat3d::rot(vec3d::Y, M_PI/2) * to_std);
            add_sub_axes(sym);
            if(nfold%2 && nfold>2) {
               sym.init(Dv, nfold/2, mat3d::rot(axis, M_PI/nfold) * to_std);
               add_sub_axes(sym);
            }
            add_sub_axes(*this);
            break;

         case D:
            dih_type = nfold%2 ? C : D;
            sym.init(dih_type, 2, mat3d::rot(vec3d::Y, M_PI/2) * to_std);
            add_sub_axes(sym);
            sym.init(dih_type, 2, mat3d::rot(axis, M_PI/nfold) *
                            mat3d::rot(vec3d::Y, M_PI/2) * to_std);
            add_sub_axes(sym);
            if(nfold%2 && nfold>2) {
               sym.init(D, nfold/2, mat3d::rot(axis, M_PI/nfold) * to_std);
               add_sub_axes(*this);
            }
            add_sub_axes(*this);
            break;

         case S:
            add_sub_axes(*this);
            break;

         case Ch:
            sym.sym_type=Cs;       // horizontal mirror
            sub_syms.insert(sym);
            add_sub_axes(*this);
            break;

         case Cv:
            if(nfold%2==0) {       // nfold even: add second vertical mirror
               sym.init(Cs, 0, mat3d::rot(axis, M_PI/nfold) *
                               mat3d::rot(vec3d::X, M_PI/2) * to_std);
               sub_syms.insert(sym);
               if(nfold>2) {
                  sym.init(Cv, nfold/2, mat3d::rot(axis, M_PI/nfold) * to_std);
                  sub_syms.insert(sym);
               }
            }
            // vertical mirror
            sym.init(Cs, 0, mat3d::rot(vec3d::X, M_PI/2) * to_std);
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
      
sch_sym_autos &sch_sym::get_autos()
{
   if(!autos.is_set())
      autos = sch_sym_autos(*this);
   return autos;
}


bool sch_sym::operator <(const sch_sym &s) const
{
   if(sym_type < s.sym_type)
      return true;
   else if(sym_type > s.sym_type)
      return false;

   if(sym_type>=C && sym_type<=S) {
      if(nfold < s.nfold)
         return true;
      else if(nfold > s.nfold)
         return false;
   }
   
   // Expensive test if same kind of symmetry group!!!
   return compare(get_trans(), s.get_trans())==-1;
}

void sch_sym_autos::init()
{
   free_vars = FREE_NONE;
   fixed_type = 0;
   fixed_trans.clear();
   for(int i=0; i<3; i++)
      rot[i] = transl[i] = 0.0;
}

sch_sym_autos::sch_sym_autos(const sch_sym &sym)
{
   init();
   const vec3d axis = vec3d::Z;
   const vec3d perp = vec3d::X;
   int type = sym.get_sym_type();
   if(type==sch_sym::unknown)    // fixed_trans empty, indicates unset;
      return;

   int nfold = sym.get_nfold();
   
   vector<vector<mat3d> > fixed;       // fixed (non-free) realignments
   fixed.push_back(vector<mat3d>(1));  // identity

   if(type==sch_sym::I || type==sch_sym::O || type==sch_sym::T ||
         type==sch_sym::D || type==sch_sym::C || type==sch_sym::C1) {
      fixed.push_back(vector<mat3d>(2));
      fixed.back()[1] = mat3d::inversion();      // reflect through origin
   }

   if(type==sch_sym::Td || type==sch_sym::T || type==sch_sym::S ||
         type==sch_sym::Ch || type==sch_sym::Cv || type==sch_sym::C) {
      fixed.push_back(vector<mat3d>(2));
      fixed.back()[1] = mat3d::rot(M_PI, 0, 0);   // flip "principal" axis
   }

   if(type==sch_sym::Dh || type==sch_sym::Dv || type==sch_sym::D ||
         type==sch_sym::Cv) {
      fixed.push_back(vector<mat3d>(2));
      fixed.back()[1] = mat3d::rot(0, 0, M_PI/nfold);// rotate base vert to edge
   }

   if((type==sch_sym::Dh || type==sch_sym::D) && nfold==2) {
      fixed.push_back(vector<mat3d>(3));
      fixed.back()[1] = mat3d::rot(vec3d(1,1,1),  2*M_PI/3); // rotate D2 axes
      fixed.back()[2] = mat3d::rot(vec3d(1,1,1), -2*M_PI/3); // rotate D2 axes
   }

   // find all combinations of transformations involving one member
   // from each set (is there an STL algorithm for this?)
   int sz = fixed.size();
   vector<unsigned int> visit_idxs(sz, 0);
   while(visit_idxs[0] < fixed[0].size()) { // until first counter rolls
      mat3d trans;
      for(int i=0; i<sz; i++)
         trans = fixed[i][visit_idxs[i]] * trans;
      fixed_trans.push_back(trans);
      //mat3d::inverse(sym.get_to_std()) * trans * sym.get_to_std() );


      // increment visit counters
      visit_idxs[sz-1]++;
      for(int j=1; j<sz; j++)
         if(visit_idxs[sz-j]>=fixed[sz-j].size()) {
            visit_idxs[sz-j] = 0;
            visit_idxs[sz-j-1]++;
         }
   }
   
   if(type==sch_sym::S || type==sch_sym::Ch)
      free_vars = FREE_ROT_PRINCIPAL;
   if(type==sch_sym::Cv)
      free_vars = FREE_TRANSL_PRINCIPAL;
   else if(type==sch_sym::C)
      free_vars = FREE_ROT_PRINCIPAL | FREE_TRANSL_PRINCIPAL; 
   else if(type==sch_sym::Cs)
      free_vars = FREE_ROT_PRINCIPAL | FREE_TRANSL_PLANE; 
   else if(type==sch_sym::Ci)
      free_vars = FREE_ROT_FULL;
   else if(type==sch_sym::C1)
      free_vars = FREE_ROT_FULL | FREE_TRANSL_SPACE; 
}

void sch_sym_autos::set_fixed(const t_set &fixed)
{
   fixed_trans.clear();
   vector<mat3d> indirect;
   fixed_trans.push_back(mat3d());
   set<mat3d>::const_iterator si;
   for(si=fixed.begin(); si!=fixed.end(); ++si) {
      iso_type iso(*si);
      if(iso.is_direct() && iso.get_rot_type()!=iso_type::rt_unit)
         fixed_trans.push_back(*si);
      else if(!iso.is_direct())
         indirect.push_back(*si);
   }
   fixed_trans.insert(fixed_trans.end(), indirect.begin(), indirect.end());
}


bool sch_sym_autos::set_fixed_type(int type, char *errmsg)
{
   int valid = false;
   if(!is_set()) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "fixed type: symmetry type not yet set");
   }
   else if(type<0 || type>=(int)fixed_trans.size()) {
      if(errmsg) {
         if(fixed_trans.size()==1)
            snprintf(errmsg, MSG_SZ, "fixed type: invalid type %d, must be 0",
                  type);
         else
            snprintf(errmsg, MSG_SZ, "fixed type: invalid type %d, must be 0-%d",
                  type, (int)fixed_trans.size()-1);
      }
   }
   else {
      fixed_type = type;
      valid = true;
   }

   return valid;
}

bool sch_sym_autos::set_rot_principal(double rot_ang, char *errmsg)
{
   int valid = free_vars & FREE_ROT_PRINCIPAL;
   if(valid) {
      rot[0] = deg2rad(rot_ang);
      rot[1] = 0.0;
      rot[2] = 0.0;
   }
   else if(errmsg)
      strcpy(errmsg, "symmetry type does not have a free principal axis rotation");
   return valid;
}

bool sch_sym_autos::set_rot_full(double rot_x, double rot_y, double rot_z,
      char *errmsg)
{
   int valid = free_vars & FREE_ROT_FULL;
   if(valid) {
      rot[0] = deg2rad(rot_x);
      rot[1] = deg2rad(rot_y);
      rot[2] = deg2rad(rot_z);
   }
   else if(errmsg)
      strcpy(errmsg, "symmetry type does not have a free full rotation");
   return valid;
}

bool sch_sym_autos::set_transl_principal(double transl0, char *errmsg)
{
   int valid = free_vars & FREE_TRANSL_PRINCIPAL;
   if(valid) {
      transl[0] = transl0;
      transl[1] = 0.0;
      transl[2] = 0.0;
   }
   else if(errmsg)
      strcpy(errmsg, "symmetry type does not have a free principal axis translation");
   return valid;
}


bool sch_sym_autos::set_transl_plane(double transl0, double transl1,
      char *errmsg)
{
   int valid = free_vars & FREE_TRANSL_PLANE;
   if(valid) {
      transl[0] = transl0;
      transl[1] = transl1;
      transl[2] = 0.0;
   }
   else if(errmsg)
      strcpy(errmsg, "symmetry type does not have a free plane translation");
   return valid;
}


bool sch_sym_autos::set_transl_space(double transl0, double transl1,
      double transl2, char *errmsg)
{
   int valid = free_vars & FREE_TRANSL_SPACE;
   if(valid) {
      transl[0] = transl0;
      transl[1] = transl1;
      transl[2] = transl2;
   }
   else if(errmsg)
      strcpy(errmsg, "symmetry type does not have a free space translation");
   return valid;
}


int sch_sym_autos::num_free_rots() const
{
   int cnt = 0;
   if(free_vars & FREE_ROT_PRINCIPAL)
       cnt = 1;
   else if(free_vars & FREE_ROT_FULL)
       cnt = 3;
   return cnt;
}
   
int sch_sym_autos::num_free_transls() const
{
   int cnt = 0;
   if(free_vars & FREE_TRANSL_PRINCIPAL)
      cnt = 1;
   else if(free_vars & FREE_TRANSL_PLANE)
      cnt = 2;
   else if (free_vars & FREE_TRANSL_SPACE)
      cnt = 3;
   return cnt;
}


mat3d sch_sym_autos::get_realignment() const
{
   mat3d trans;
   if(is_set())
      trans = fixed_trans[fixed_type];
   
   if(free_vars & FREE_ROT_PRINCIPAL)
      trans = mat3d::rot(vec3d::Z, rot[0]) * trans;
   else if(free_vars & FREE_ROT_FULL)
      trans = mat3d::rot(rot[0], rot[1], rot[2]) * trans;

   if(free_vars & FREE_TRANSL_PRINCIPAL)      // z-axis
      trans = mat3d::transl(vec3d(0, 0, transl[0])) * trans;
   else if(free_vars & FREE_TRANSL_PLANE)     // xy-plane
      trans = mat3d::transl(vec3d(transl[0], transl[1], 0)) * trans;
   else if (free_vars & FREE_TRANSL_SPACE)
      trans = mat3d::transl(vec3d(transl[0], transl[1], transl[2])) * trans;

   return trans;
}
      
bool sch_sym_autos::set_realignment(const char *realign, char *errmsg)
{
   char errmsg2[MSG_SZ];
   vector<double> vars;
   
   char *realign_cpy = copy_str(realign);
   char *p = strchr(realign_cpy, ':');
   if(p)
      *p = '\0';                 // terminator at first ':'
   
   bool valid = true;
   int type=0;
   if(*realign=='\0' || (p && p == realign_cpy)) // empty, so set default of 0
      type = 0;
   else {
      if(!read_int(realign_cpy, &type, errmsg2)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "fixed type: %s", errmsg2);
         valid =  false;
      }
   }

   if(valid && !set_fixed_type(type, errmsg))
      valid = false;

   if(valid && p && *(p+1)) {    // colon found and characters afterwards
      if(!read_double_list(p+1, vars, errmsg2, 0, ":")) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "free variables: %s", errmsg2);
         valid = false;
      }
   }

   free(realign_cpy);
   if(!valid)
      return false;

   int rot_cnt = num_free_rots();         // number of rotation angles 
   int transl_cnt = num_free_transls();   // number of translation distances 

   int vars_sz = vars.size();
   if(vars_sz && vars_sz!=rot_cnt && vars_sz!=rot_cnt+transl_cnt) {
      if(errmsg) {
         string msg = msg_str("free variable list: %d numbers given",
               vars.size());
         if(rot_cnt+transl_cnt==0)
            msg += " but there are no free variables";
         else if(rot_cnt>0 && transl_cnt>0)
            msg += msg_str(", must give %d (rotation) or %d (rotation then "
                  "translation) colon separated numbers",
                  rot_cnt, rot_cnt+transl_cnt);
         else if(rot_cnt>0)
            msg += msg_str(", must give %d (rotation) colon separated numbers",
                  rot_cnt);
         else //(transl_cnt>0)
            msg += msg_str(", must give %d (translation) colon separated "
                  "numbers", rot_cnt);

         strncpy(errmsg, msg.c_str(), MSG_SZ);
      }
      return false;
   }

   if(vars_sz) {
      if(free_vars & FREE_ROT_PRINCIPAL)
         set_rot_principal(vars[0]);
      else if(free_vars & FREE_ROT_FULL)
         set_rot_full(vars[0], vars[1], vars[2]);
   }

   if(vars_sz>rot_cnt) {
      if(free_vars & FREE_TRANSL_PRINCIPAL)
         set_transl_principal(vars[rot_cnt+0]);
      else if(free_vars & FREE_TRANSL_PLANE)
         set_transl_plane(vars[rot_cnt+0], vars[rot_cnt+1]);
      else if (free_vars & FREE_TRANSL_SPACE)
         set_transl_space(vars[rot_cnt+0], vars[rot_cnt+1], vars[rot_cnt+2]);
   }

   return true;
}

