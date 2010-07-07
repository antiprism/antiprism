/*
   Copyright (c) 2003-2008, Adrian Rossiter

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

/*!\file polygons.h
   \brief Generate polyhedra based on polygons.
*/

#ifndef POLYGONS_H
#define POLYGONS_H


#include <string.h>

#include "geom.h"
#include "transforms.h"
#include "math_utils.h"

///Make a uniform %polygon
class polygon {
   protected:
      int num_sides; ///< The number of sides to the %polygon.
      int fraction;  ///< The value m in the {n/m} description.
      int parts;     ///< The number of parts (the %polygon may be a compound).
      double radius; ///< The circumradius.
      double height; ///< The height, of a %polygon based polyhedron.
      
   public:
      ///Constructor
      /**\param sides the number of sides to the %polygon.
       * \param fract the value m in the {n/m} description. */
      polygon(int sides, int fract=1);
      
      ///Destructor
      virtual ~polygon() {}
      
      ///Set the circumradius
      /**\param r the circumradius. */
      void set_radius(double r) { radius = r; }

      ///Set the inradius
      /**\param r the inradius. */
      void set_inradius(double r) { radius = r / cos(angle()/2);}

      ///Set the edge (%polygon side) length
      /**\param len the edge length. */
      void set_edge(double len) { radius = 0.5*len/sin(angle()/2); }

      ///Get the edge (polygon side) length
      /**\return len the edge length. */
      double get_edge() { return 2*radius*sin(angle()/2); }

      ///Get the angle that an edge makes at the centre
      /**\return the angle, in radians. */
      double angle() { return 2*M_PI*fraction/num_sides; }

      ///Get the number of sides of the %polygon.
      /**\return the number of sides. */
      int get_num_sides() { return num_sides; }

      ///Add a %polygon aligned with the xy plane and at a given z-height
      /**\param geom the geometry to add the %polygon to.
       * \param ht the height on the z-axis to place the %polygon. */
      void add_polygon(geom_if &geom, double ht);

      ///Make a %polygon based polyhedron.
      /**\param geom a geometry to return the polyhedron. */
      void make_poly(geom_if &geom);
     
      ///Set the height
      /**\param ht the height.
       * \param msg a string with length at least \c MSG_SZ to hold
       * the error message if the height was not valid.
       * \return \c true if the height was valid, otherwise \c false and
       * \c msg contains the error messge. */
      virtual bool set_height(double ht, char *msg=0)
         { height = ht; msg=0; return true;}

      ///Set the edge length of the non-polygon edges.
      /**These are the vertical edges of a prism, the slanting edges
       * of a pyramid, etc.
       * \param len2 the edge length of the non-polygon edges.
       * \param msg a string with length at least \c MSG_SZ to hold
       * the error message if the edge length was not valid.
       * \return \c true if the edge length was valid, otherwise \c false and
       * \c msg contains the error messge. */
      virtual bool set_edge2(double len2, char *msg=0)
         {  len2=0;
            if(msg) strcpy(msg, "cannot set slanting edge for this polyhedron");
            return false; }
      
      ///Set the sub-type of the %polygon based polyhedron.
      /**Some %polygon based polyhedra come in several forms, and setting
       * the sub-type can select a particular form.
       * \param typ the sub-type number.
       * \param msg a string with length at least \c MSG_SZ to hold
       * the error message if the edge length was not valid.
       * \return \c true if the sub-type was valid, otherwise \c false and
       * \c msg contains the error messge. */
      virtual bool set_subtype(int typ, char *msg=0)
         {  typ=0;
            if(msg) strcpy(msg, "this type of polyhedron has no sub-types");
            return false; }

      ///Make a part of (or a complete) polygon-based polyhedron
      /**Make a non-compound polyhedron, using \c num_sides and
       * \c fraction for {n/m}. If \c parts is greater than \c 1
       * then polygon::make_poly will make a compound by repeating this
       * polyhedron \c parts times.
       * \param geom a geometry to return the polyhedron. */
      virtual void make_poly_part(geom_if &geom) { add_polygon(geom, height); }
};


///Make a %dihedron
class dihedron: public polygon {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      dihedron(int sides, int fract=1) : polygon(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      dihedron(polygon &pgon) : polygon(pgon) { height = radius;}
      
      bool set_edge2(double e2, char * /*msg*/ =0) { height = e2; return true; }
      void make_poly_part(geom_if &geom);
}; 


///Make a %prism
class prism: public polygon {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      prism(int sides, int fract=1) : polygon(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      prism(polygon &pgon) : polygon(pgon) { height = radius;}

      bool set_edge2(double e2, char * /*msg*/ =0) { height = e2; return true; }
      void make_poly_part(geom_if &geom);
}; 


///Make an %antiprism
class antiprism: public polygon {
   private:
      double twist_angle;
      bool output_trapezohedron;
      void make_trapezo_part(geom_if &geom);
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      antiprism(int sides, int fract=1) : polygon(sides, fract),
         twist_angle(0.0), output_trapezohedron(false) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      antiprism(polygon &pgon) : polygon(pgon),
         twist_angle(0.0), output_trapezohedron(false)
         { height = radius;}

      bool set_edge2(double e2, char *msg=0);
      void set_output_trapezohedron(bool flag=true)
         { output_trapezohedron=flag; }
      void set_twist_angle(double ang) { twist_angle=ang; }
      void make_poly_part(geom_if &geom);
}; 


///Make a %snub-antiprism
class snub_antiprism: public polygon {
   private:
      int type;
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      snub_antiprism(int sides, int fract=1) : polygon(sides, fract), type(0) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      snub_antiprism(polygon &pgon) : polygon(pgon), type(0) { height = radius;}
     
      bool set_edge2(double e2, char *msg=0);
      bool set_height(double /*h*/, char *msg=0);
      bool set_subtype(int typ=0, char *msg=0);
      void make_poly_part(geom_if &geom);
}; 


///Make a %pyramid
class pyramid: public polygon {
   private:
      bool output_trapezohedron;

   public:
      ///Constructor
      /**\param sides the number of sides to the base %polygon.
       * \param fract the value m in the {n/m} description. */
      pyramid(int sides, int fract=1) : polygon(sides, fract),
         output_trapezohedron(false) {}
      
      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      pyramid(polygon &pgon) : polygon(pgon), output_trapezohedron(false)
         { height = radius;}
      
      bool set_edge2(double e2, char *msg=0);
      void set_output_trapezohedron(bool flag=true)
         { output_trapezohedron=flag; }
      void make_poly_part(geom_if &geom);
}; 


///Make a %dipyramid
class dipyramid: public pyramid {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      dipyramid(int sides, int fract=1) : pyramid(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      dipyramid(polygon &pgon) : pyramid(pgon) {}
      
      void make_poly_part(geom_if &geom);
}; 


///Make a %cupola
class cupola: public polygon {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
     cupola(int sides, int fract=1) : polygon(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
     cupola(polygon &pgon) : polygon(pgon) { height = radius;}
     
     bool set_edge2(double e2, char *msg=0);
     void make_poly_part(geom_if &geom);
}; 


///Make an %orthibicupola
class orthobicupola: public cupola {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      orthobicupola(int sides, int fract=1) : cupola(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      orthobicupola(polygon &pgon) : cupola(pgon) {}

      void make_poly_part(geom_if &geom);
};


///Make a %gyrobicupola
class gyrobicupola: public cupola {
   public:
      ///Constructor
      /**\param sides the number of sides to the base-polygon.
       * \param fract the value m in the {n/m} description. */
      gyrobicupola(int sides, int fract=1) : cupola(sides, fract) {}

      ///Constructor
      /**\param pgon %polygon to base the polyhedron on. */
      gyrobicupola(polygon &pgon) : cupola(pgon) {}

      void make_poly_part(geom_if &geom);
};

///Make a uniform model of a polygon-based polyhedron.
/**The model has all its edges set to one and its faces are regular. For
 * prism's, antiprisms and pyramids the resulting model wil be uniform.
 * \param geom a geometry to return the model.
 * \param pgon a polygon-derived object of the polyhedron type required;
 * \return true if the edge types could be set to length \c 1.0,
 * otherwise \c false. */
template <class T> bool uni_pgon(geom_if &geom, T pgon)
{ 
   pgon.set_edge(1.0);
   bool ret = pgon.set_edge2(1.0, 0);
   pgon.make_poly(geom);
   return ret;
}

int make_resource_pgon(geom_if &geom, string pname, char *errmsg=0);
      

// Inline function definitions

inline polygon::polygon(int sides, int fract) :
   fraction(fract), radius(1.0), height(0)
{
   parts = gcd(sides, fract);
   num_sides = sides/parts;
   fraction = fract/parts;
}


#endif // POLYGONS_H

