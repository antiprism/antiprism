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
      
      ///Get the underlying transformation set
      /** return The underlying transformation set */
      const set<mat3d> &get_trans() const { return trans; }
      
      ///Get the underlying transformation set
      /** return The underlying transformation set */
      set<mat3d> &get_trans() { return trans; }
      
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

///Compare two t_sets for order
/**Returns -1, 0, or 1 to indicate less, equal or greater. Order by
 * size, and if equal compare matrices sequentially until not equal,
 * and if all equal return 0.
 * \param t0 first t_set
 * \param t1 second t_set
 * \return  -1, 0, or 1 to indicate t0 is less, equal or greater than t1 */
int compare(const t_set &t0, const t_set &t1);


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

      ///Check if isometry is direct.
      /** \return \c true if the isometry is direct, otherwise \cfalse */
      bool is_direct() const { return rot_type>rt_none && rot_type<rt_inv; }

      ///Dump
      /**Print the object data to \c stdout for debugging. */
      void dump() const;
      
};


///Get sets of elements that are equivalent under a set of transformations
/**\param geom the geometry.
 * \param ts the set of transfromations to apply.
 * \param equiv_sets vectors of sets of equivalent elements for
 * vertices (0), edges (1) and faces (2). */
void get_equiv_elems(const geom_if &geom, const t_set &ts,
      vector<vector<set<int> > > *equiv_sets);

class sch_axis;
class sch_sym;

///Class for a set of automorphisms: transformations maintaining symmetry alignment 
class sch_sym_autos
{
   private:
      vector<mat3d> fixed_trans;
      unsigned int free_vars;
      int fixed_type;
      double rot[3];
      double transl[3];

      void init();

   public:
      enum {FREE_NONE=0, FREE_ROT_PRINCIPAL=1, FREE_ROT_FULL=2,
           FREE_TRANSL_PRINCIPAL=4, FREE_TRANSL_PLANE=8, FREE_TRANSL_SPACE=16 };

      ///Constructer
      sch_sym_autos() { init(); }

      ///Constructer
      /**\param sym symmetry group to transform */
      sch_sym_autos(const sch_sym &sym);

      ///Check if has been set
      bool is_set() const { return fixed_trans.size(); }

      ///Set fixed transformation type
      /**\param type the index number of the type (0 always the identity)
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the fixed transformation index is valid,
       * otherwise false and the error is detailed in \a errmsg. */
      bool set_fixed_type(int type, char *errmsg=0);

      ///Set rotation about principal axis
      /**\param rot_ang rotation angle about prinicipal axis, in degrees
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the free variable could be set, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_rot_principal(double rot_ang, char *errmsg=0);

      ///Set rotation about origin
      /**\param rot_x rotation about x-axis, in degrees
       * \param rot_y rotation about y-axis, in degrees
       * \param rot_z rotation about z-axis, in degrees
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the free variable could be set, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_rot_full(double rot_x, double rot_y, double rot_z,
            char *errmsg=0);

      ///Set translation distance along principal direction
      /**\param transl0 translation distance
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the free variable could be set, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_transl_principal(double transl0, char *errmsg=0);

      ///Set translation distances along directions which span a plane
      /**\param transl0 first translation distance
       * \param transl1 second translation distance
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the free variable could be set, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_transl_plane(double transl0, double transl1, char *errmsg=0);

      ///Set translation distances along directions which span all space
      /**\param transl0 first translation distance
       * \param transl1 second translation distance
       * \param transl2 third translation distance
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the free variable could be set, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_transl_space(double transl0, double transl1, double transl2,
            char *errmsg=0);

      /// Get  fixed realignment transformations
      /**\return The fixed realigment transformations.*/
      const vector<mat3d> &get_fixed() const { return fixed_trans; }

      /// Set  fixed realignment transformations
      /**\param fixed the fixed realigment transformations.*/
      void set_fixed(const t_set &fixed);

      /// Number of free rotation variables
      /**\return Number of free rotation variables (0, 1 or 3).*/
      int num_free_rots() const;
      
      /// Number of free translation variables
      /**\return Number of free translation variables (0 to 3).*/
      int num_free_transls() const;
      
      ///Get realignment transformtion from current settings
      /**\return Realignment transformation */
      mat3d get_realignment() const;
      
      ///Set and get realignment transformation from a string of numbers
      /**\param realign fixed type number and free variable list as
       * colon separated string
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the realignment was valid, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_realignment(const char *realign, char *errmsg=0);
};


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
      mutable set<sch_sym> sub_syms;
      mutable set<sch_axis> axes;
      mutable set<vec3d> mirrors;
      sch_sym_autos autos;
      mat3d to_std;

      void add_sub_axes(const sch_sym &sub) const;
      void find_full_sym_type(const set<sch_axis> &full_sym);

      ///Set the symmetry type for the axis
      /**\param type the symmetry type of the axis as the Schoenflies
       * identifier. */
      void set_sym_type(int type) { sym_type=type; }
      
      ///Set the n-fold order of the axis.
      /**If the symmetry type is \c sym_S then the axis has
       * rotational n/2-fold symmetry.
       * \param n the n-fold order of the axis. */
      void set_nfold(int n) { nfold = n; }

      ///Set the tranformation to standard symmetry type.
      /**\param trans the transformation that carries an object with the
       * symmetry type onto the standard set of symmetries for that type. */
      void set_to_std(const mat3d &trans) { to_std=trans; }



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
      sch_sym(const string &name, const mat3d &pos=mat3d(), char *errmsg=0);

       ///Constructor
      /**\param sym_ax symmetry type in symmetry axis form.
       * \param cent coordinates of symmetry centre */
      sch_sym(const sch_axis &sym_ax, const vec3d &cent=vec3d::zero);

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
      bool init(const string &name, const mat3d &pos=mat3d(), char *errmsg=0);

       ///Initialiser
      /**\param sym_axis symmetry type in symmetry axis form.
       * \param cent coordinates of symmetry centre */
      void init(const sch_axis &sym_axis, const vec3d &cent=vec3d::zero);

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

      ///Get the symmetry subgroups
      /**Only one example is included from each conjugacy class
       * \return The symmetry subgroups. */
      const set<sch_sym> &get_sub_syms() const;

      ///Get the Euclidean outer automorphisms
      /**These map the symmetry transformations onto themselves, with one
       * example included from each conjugacy class, the first value 
       * in the returned transformations is always the identity.
       * \return The symmetry automorphisms. */
      sch_sym_autos &get_autos();

      ///Get a symmetry subgroup
      /**\param sub_sym the symmetry subgroup
       * \param conj_type use to select from inequivalent (non-conjugate)
       * subgroups.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return the sub-symmetry, if the symmetry type is 'unknown' the
       * error is detailed in \a errmsg. */
      sch_sym get_sub_sym(const sch_sym &sub_sym, int conj_type=0,
            char *errmsg=0) const;
      
      
      ///Get a symmetry subgroup
      /**\param sub the symmetry subgroup name and conjugation number,
       * separated by a comma
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return the sub-symmetry, if the symmetry type is 'unknown' the
       * error is detailed in \a errmsg. */
      sch_sym get_sub_sym(const string &sub, char *errmsg) const;


      ///Get the axes or mirrors.
      /**\return The symmetry axes or mirrors. */
      const set<sch_axis> &get_axes() const;

      ///Get the mirrors.
      /**\return The mirror normals. */
      const set<vec3d> &get_mirrors() const;

      ///Is an inversion a symmetry.
      /**\return \c true if an inversion is a symmetry, otherwise \c false */
      bool has_inversion_symmetry() const;

      ///Check if a valid symmetry type is set.
      /**\return \c true if set to valid type, otherwise \c false */
      bool is_set() const { return sym_type!=unknown; }

      // Less than.
      /* Only for containers that need it
       * \param s the symmetry axis to compare for less than.
       * \return \c true if this axis is less than \arg s, otherwise \c false.*/
      bool operator <(const sch_sym &s) const;


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


