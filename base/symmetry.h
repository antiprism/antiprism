/*
   Copyright (c) 2003-2008, Adrian Rossiter

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


/*!\file symmetry.h
 *\brief symmetry and transformation handling
*/


#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <vector>
#include <set>
#include <string>

#include <math.h>
#include "mat3d.h"

using std::vector;
using std::set;
using std::string;

/// A set of transformations
class t_set
{
   private:
      set <mat3d> trans;

   public:
      ///A constant iterator over the set
      typedef set<mat3d>::const_iterator const_iterator;
      
      ///Add a transformation
      /**\param tr transformation to add */
      t_set &add(const mat3d &tr) { trans.insert(tr); return *this; }
      
      ///Clear all transformations
      /**Remove all the transformations */
      t_set &clear() { trans.clear(); return *this; }

      ///Set to direct product of two sets of transformations
      /**Each transformation of one set is combined with every
       * transformation from the other set.
       * \param s1 one of the sets.
       * \param s2 the other set.
       * \return A reference to this object with the direct product
       * of transformations set. */
      t_set &product(const t_set& s1, const t_set &s2);

      ///Set to direct product of current transformations with another set of transformations
      /**Each transformation of current set is combined with every
       * transformation from the other set.
       * \param s the other set.
       * \return A reference to this object with the direct product
       * of transformations set. */
      t_set &product_with(const t_set& s);

      ///Set transformations to their conjugates by a given transformation.
      /**If the given transformation is M then each transformation
       * T is set to MTM^-1.
       * \param m the transformation to form the conjugates with.
       * \return A reference to this object with the conjugates set. */
      t_set &conjugate(const mat3d &m);

      ///Intersection of two sets of transformations.
      /**Find the transformations that are common to both sets of
       * transformations.
       * \param s1 one of the sets.
       * \param s2 the other set.
       * \return A reference to this object with common transformations set. */
      t_set &intersection(const t_set &s1, const t_set &s2);

      ///Subtract from current transformations any common to a second set of transformations.
      /**After this operation the two sets will have no transformations
       * in common.
       * \param s the second set.
       * \return A reference to this object with the common transformations
       * removed. */
      t_set &subtract(const t_set &s);

      ///Form the left coset with a given transformation.
      /**If the given transformation is M then each transformation
       * T is set to MT in the coset.
       * \param m the transformation to form the left coset with.
       * \return The coset with \arg m. */
      t_set lcoset(const mat3d &m) const;

      ///The minimum set of transformations that will generate one group of transformations from another
      /**\param tr_whole the target group of transformations.
       * \param tr_part the group of transformations that will be used to
       * generate the target group.
       * \param pos A transformation that will align \arg tr_part with
       * \arg tr_whole.
       * \return A reference to this object containing a minimum set of
       * transformations that when applied to a geometry with symmetry
       * \arg tr_part will create a compound with symmetry \arg tr_whole. */
      t_set &min_set(const t_set &tr_whole, const t_set &tr_part,
            const mat3d &pos=mat3d());

      ///Begin, for iterating through the set.
      /**\return A constant iterator pointing to the first transformation. */
      const_iterator begin() const { return trans.begin(); }

      ///End, for iterating through the set.
      /**\return A constant iterator pointing to one past the last
       * transformation. */
      const_iterator end() const { return trans.end(); }

      ///Size, the number of transformations.
      /**\return The number of transformations. */
      size_t size() const { return trans.size(); }
      
      //const set<mat3d> &get_trans() const { return trans; }
      
      ///Check if there are transformations in the set.
      /**\return \c true if set, otherwise \c false */
      bool is_set() const { return trans.size()>0; }

};

///Form direct product of two sets of transformations
/**Each transformation from one set is combined with every transformation
 * from the other set.
 * \param s1 one of the sets.
 * \param s2 the other set.
 * \return The direct productof the two sets. */
t_set operator *(const t_set &s1, const t_set &s2);

///Form direct product of two sets of transformations
/**Each transformation of the first set is combined with every
 * transformation from the second set and the result added to the
 * first set.
 * \param s1 the first set.
 * \param s2 the second set.
 * \return A reference to s1 with the direct product of transformations set. */
t_set &operator *=(t_set &s1, const t_set &s2);

///Add a transformation
/**\param s a set of transformations.
 * \param m the transformation to add.
 * \return A copy of \arg s with \arg m added. */
t_set operator +(const t_set &s, const mat3d &m);

///Add a transformation
/**\param s a set of transformations.
 * \param m the transformation to add.
 * \return A refernce to \arg s with \arg m added. */
t_set &operator +=(t_set &s, const mat3d &m);


///Class for an isometry
class iso_type
{
   public:
      ///Rotation type
      enum {rt_none=0, rt_unit, rt_rot, rt_inv, rt_refl, rt_rot_refl};
   private:
      int rot_type;
      vec3d axis;
      double ang;
      vec3d transl;

      bool is_isometry(mat3d m) const;
      void normalise();
      
   public:
      ///Constructor
      iso_type() : rot_type(rt_none) {}
      
      ///Constructor
      /**Set up with the details for a particular transformation matrix.
       * \param m the transformation matrix. */
      iso_type(const mat3d &m) { init(m); }

      ///Initialise with a transformation matrix.
      /**Set up with the details for a particular transformation matrix.
       * \param m the transformation matrix.
       * \return A reference to this symmetry axis. */
      iso_type &init(mat3d m);

      ///Get the rotation type
      /**A return value of \c rt_none indicates that the object has not
       * been initialised with an isometry.
       * \return The rotation type. */
      int get_rot_type() const { return rot_type; }

      ///Get principal axis.
      /**The axis will be unset if the rotation type doesn't have an axis.
       * \return The axis. */
      vec3d get_axis() const { return axis; }
      
      ///Get rotation angle.
      /**The angle will be 0 if the rotation type doesn't have an angle.
       * \return The angle. */
      double get_ang() const { return ang; }
      
      ///Get translation.
      /**The translation will be unset if the object has not been
       * initialised with an isometry.
       * \return The translation. */
      vec3d get_transl() const { return transl; }

      ///Dump
      /**Print the object data to \c stdout for debugging. */
      void dump() const;
      
};

class sch_axis;

///Class for symmetry elements in Schoenflies notation.
class sch_sym
{
   public:
      ///Schoenflies symmetry type identifiers
      enum { unknown, C1, Ci, Cs, C, Cv, Ch, D, Dv, Dh, S,
         T, Td, Th, O, Oh, I, Ih };

   private:
      int sym_type;
      int nfold;
      mutable set<sch_axis> axes;
      mutable set<vec3d> mirrors;
      mat3d to_std;

      void find_full_sym_type(const set<sch_axis> &full_sym);

   public:
      ///Constructor
      /**Find the symmetry in Schoenflies notation.
       * \param geom geometry to find the symmetry type for.
       * \param equiv_sets for the vertices, edges and faces set up a
       * vector of sets of equivalent elements. */
      sch_sym(const geom_if &geom, vector<vector<set<int> > > *equiv_sets=0);

      ///Constructor
      /**Find the symmetry in Schoenflies notation.
       * \param ts group of transformations to find the symmetry type for.*/
      sch_sym(const t_set &ts);

       ///Constructor
      /**\param type the symmetry type.
       * \param n n-fold order, for a principal axis, otherwise 0.
       * \param pos transformation that carries an object of this
       * symmetry type onto the standard symetry type.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message. */
      sch_sym(int type=0, int n=0, const mat3d &pos=mat3d(), char *errmsg=0);

       ///Constructor
      /**\param name name of symmetry type to set up.
       * \param pos transformation that carries an object of this
       * symmetry type onto the standard symetry type.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message. */
      sch_sym(string name, const mat3d &pos=mat3d(), char *errmsg=0);

      ///Initialiser
      /**Find the symmetry in Schoenflies notation.
       * \param geom geometry to find the symmetry type for.
       * \param equiv_sets for the vertices, edges and faces set up a
       * vector of sets of equivalent elements
       * \return true if the symmetry type could be determined,
          * otherwise false.*/
      bool init(const geom_if &geom, vector<vector<set<int> > > *equiv_sets=0);

       ///Initialiser
      /**\param type the symmetry type.
       * \param n n-fold order, for a principal axis, otherwise 0.
       * \param pos transformation that carries an object of this
       * symmetry type onto the standard symetry type.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the symmetry specification was valid, otherwise false
       * and the error is detailed in \a errmsg. */
      bool init(int type, int n=0, const mat3d &pos=mat3d(), char *errmsg=0);

       ///Initialiser
      /**\param name name of symmetry type to set up.
       * \param pos transformation that carries an object of this
       * symmetry type onto the standard symetry type.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the symmetry type name was valid, otherwise false
       * and the error is detailed in \a errmsg. */
      bool init(string name, const mat3d &pos=mat3d(), char *errmsg=0);

      ///Get the symmetry type
      /**\return the symmetry type of the axis as the Schoenflies
       * identifier. */
      int get_sym_type() const { return sym_type; }

      ///Get the n-fold order of the axis with highest n.
      /**This value is used in conjunction with \c get_sym_type()
       * to describe principal axis symmetry types. If the symmetry
       * type is \c sym_S and the axis has rotational n-fold symmetry
       * then 2n is returned.
       * \return the n-fold order of the axis. */
      int get_nfold() const { return nfold; }
      
      ///Get the Schoenflies symbol for the symmetry type
      /**\return The symbol. */
      string get_symbol() const;

      ///Get the tranformation to standard symmetry type.
      /**\return The transformation that carries an object with the symmetry
       * type onto the standard set of symmetries for that type. */
      mat3d get_to_std() const { return to_std; }

      ///Get the symmetry transformations
      /**\param ts to return the set of symmetry transformations for this
       * symmetry type. */
      t_set &get_trans(t_set &ts) const;

      ///Get the symmetry transformations
      /**\return The set of symmetry transformations for this symmetry type. */
      t_set get_trans() const;

      ///Get the axes or mirrors.
      /**\return The symetry axes or mirrors. */
      const set<sch_axis> &get_axes() const;

      ///Get the mirrors.
      /**\return The mirror normals. */
      const set<vec3d> &get_mirrors() const;

      ///Check if a valid symmetry type is set.
      /**\return \c true if set to valid type, otherwise \c false */
      bool is_set() const { return sym_type!=unknown; }

};

///Class for a principal axis or mirror symmetry in Schoenflies notation.
class sch_axis
{
   private:
      int sym_type;
      vec3d axis;
      vec3d perp;
      long nfold;
      
   public:
      ///Constructor
      /**Set up with the smallest symmetry type that includes a
       * particular transfomation.
       *\param m the transformation. */
      sch_axis(const mat3d &m = mat3d());

      ///Less than.
      /**Only for containers that need it
       * \param s the symmetry axis to compare for less than.
       * \return \c true if this axis is less than \arg s, otherwise \c false.*/
      bool operator <(const sch_axis &s) const;

      ///Set the principal axis
      /**\param ax the axis. */
      void set_axis(vec3d ax) { axis = ax; }

      ///Get the principal axis
      /**\return the axis. */
      const vec3d &get_axis() const { return axis; }

      ///Set the perpendicular axis
      /**\param perp_ax the axis. */
      void set_perp(vec3d perp_ax) { perp = perp_ax; }

      ///Get the perpendicular axis
      /**Not all symmetry types have a perpendicular axis.
       * \return The perpendicular axis. */
      const vec3d &get_perp() const { return perp; }

      ///Set the symmetry type for the axis
      /**\param type the symmetry type of the axis as the Schoenflies
       * identifier. */
      void set_sym_type(int type) { sym_type=type; }
      
      ///Get the symmetry type for the axis
      /**\return the symmetry type of the axis as the Schoenflies
       * identifier. */
      int get_sym_type() const { return sym_type; }

      ///Check for horizontal mirror
      /**\return \c true if there is a horizontal mirror,
       * otherwise \c false. */
      bool has_horz_mirror() const { return ( sym_type==sch_sym::Cs || 
            sym_type==sch_sym::Ch || sym_type==sch_sym::Dh); }


      ///Set the n-fold order of the axis.
      /**If the symmetry type is \c sym_S then the axis has
       * rotational n/2-fold symmetry.
       * \param n the n-fold order of the axis. */
      void set_nfold(int n) { nfold = n; }

      ///Get the n-fold order of the axis.
      /**If the symmetry type is \c sym_S and the axis has
       * rotational n-fold symmetry then 2n is returned.
       * \return the n-fold order of the axis. */
      int get_nfold() const { return nfold; }
      
      ///Dump
      /**Print the object data to \c stdout for debugging. */
      void dump() const;
};




#endif // SYMMETRY_H


