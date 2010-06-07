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

/*!\file col_map.h
   \brief A colour map class
*/

#ifndef COL_MAP_H
#define COL_MAP_H

#include <limits.h>
#include <string>
#include <map>

using std::string;
using std::map;

///A colour map class
class color_map
{
   private:
      int shift;    // the first entry starts this many places into the map.
      int step;     // step this many places with each index increment.
      int wrap;     // wrap index numbers back to zero at this index size.

   protected:
      ///Copy parameters
      /**\param cmap the colour map to copy parameters from. */
      void copy_params(const color_map &cmap);

   public:
      const static int max_map_sz = INT_MAX;  // The largest size of a map  

      ///Constructor
      color_map(): shift(0), step(1), wrap(0) {}

      ///Destructor
      virtual ~color_map() {}
      
      ///Initialise from a string
      /** \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *map_name, char *errmsg =0);

      ///Initialise from a string, stripping parameters from the end
      /** \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool init_strip(char *map_name, char *errmsg=0);

      ///Get a copy of the map with any with the values indexed by order
      /** \return a pointer to the dynamically allocated condensed copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the copy failed. */
      virtual color_map *get_condensed() { return clone(); }

      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      virtual color_map *clone() { return new color_map(*this); }

      ///The effective size of the map
      /**The effective size of a map is one greater than the highest
       * index number in the map. It is the size of the smallest
       * map (sequential, starting at 0) that will include all
       * the entries of the map.
       * \return The effective size */
      virtual int effective_size() const { return max_map_sz; };

      ///Get the map shift
      /**Lookup of colour values is determined by (shift + index*step)%wrap
       * \return the shift */
      int get_shift() const { return shift; }

      ///Set the map shift
      /**Lookup of colour values is determined by (shift + index*step)%wrap
       * \param shft the shift */
      void set_shift(int shft) { shift=shft; }

      ///Get the map step
      /**Look up of colour values is determined by (shift + index*step)%wrap
       * \return the step */
      int get_step() const { return step; }

      ///Set the map shift
      /**Look up of colour values is determined by (shift + index*step)%wrap
       * \param stp the step */
      void set_step(int stp) { step=stp; }

      ///Get the map wrap size
      /**Look up of colour values is determined by (shift + index*step)%wrap
       * a value of 0 indicates no wrapping.
       * \return the wrap */
      int get_wrap() const { return wrap; }

      ///Set the map wrap size
      /**Look up of colour values is determined by (shift + index*step)%wrap
       * a value of 0 indicates no wrapping, a negative value indicates
       * wrapping after the highest index number.
       * \param wrp the wrap */
      void set_wrap(int wrp=-1);

      ///Get effective index
      /**The effective index is (shift + index*step)%wrap
       * \param idx the index
       * \return The effective index. */
      int get_effective_index(int idx) const;

      ///Get the colour value for an index number.
      /**
       * \return The colour. */
      virtual col_val get_col(int /*idx*/) { return col_val(); }

      ///Cycle the map colours
      /**Each colour index number is mapped to the previous colour
       * value in the colour map. */
      void cycle_colors();

};


///A colour map that remaps index numbers 
class color_map_remap: public color_map
{
   public:
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_remap(*this); }

      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx) { return get_effective_index(idx); }
};


///A colour map using a range 
class color_map_range: public color_map
{
   private:
      int map_sz;

   protected:
      vector<double> ranges[4];
      void (col_val::*set_func)(double, double, double, double);

   public:
      ///Initialise from a string
      /** \param name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *name, char *errmsg=0);

      ///Get a copy of the map
      /**\return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range(*this); }

      ///Set a range
      /**\param idx the index number of the component (0-3 for RGBA or HSVA)
       * \param range the range to set.
       * \return true if the range was valid, else false and the
       * range was not changed. */
      virtual bool set_range(int idx, vector<double> range);

      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx);

      ///Set the map size
      /**\param sz the number of entries in the map. */
      void set_map_sz(int sz) { map_sz = sz; }

      ///Get the map size
      /**\return the number of entries in the map. */
      int get_map_sz() { return map_sz; }
      
      ///The effective size of the map
      /**The effective size of a map is one greater than the highest
       * index number in the map. It is the size of the smallest
       * map (sequential, starting at 0) that will include all
       * the entries of the map.
       * \return The effective size */
      virtual int effective_size() const { return map_sz; };

};

///A colour map using values in an HSVA range 
class color_map_range_hsv: public color_map_range
{
   public:
      ///Initialise from a string
      /** \param name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *name, char *errmsg=0);
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range_hsv(*this); }

};


///A colour map using values in an RGBA range 
class color_map_range_rgb: public color_map_range
{
   public:
      ///Destructor
      ~color_map_range_rgb() {}

      ///Initialise from a string
      /** \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *map_name, char *errmsg=0);
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range_rgb(*this); }

};

///A colour map using random values in a range 
class color_map_range_rand: public color_map_range
{
   public:
      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx);
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range_rand(*this); }

};


///A colour map using random values in an HSVA range 
class color_map_range_rand_hsv: public color_map_range_rand
{
   public:
      ///Initialise from a string
      /** \param name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message. */
      color_map_range_rand_hsv(const char *name="", char *errmsg=0)
         { color_map_range_rand_hsv::init(name, errmsg); }

      ///Initialise from a string
      /** \param name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *name, char *errmsg=0);
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range_rand_hsv(*this); }

};


///A colour map using random values in an RGBA range 
class color_map_range_rand_rgb: public color_map_range_rand
{
   public:
      ///Initialise from a string
      /** \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *map_name, char *errmsg=0);
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_range_rand_rgb(*this); }

};

///A colour map with a good spread of colours 
class color_map_spread: public color_map_range
{
   public:
      ///Initialise from a string
      /** \param name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *name, char *errmsg=0);

      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_spread(*this); }

      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx);
};



///A colour map with the mappings held in a map
class color_map_map : public color_map
{
   private:
      map<int, col_val> cmap;              // color map

      bool generate_map(int map_sz, const vector<double> &comp0,
            const vector<double> &comp1, const vector<double> &comp2,
            const vector<double> &comp3,
            void (col_val::*set_func)(double, double, double, double));
   public:
      // Constructor
      color_map_map() : color_map() {}

      ///Initialise from a file
      /**The colour map can be in the Antiprism, GIMP or Fractint format. If
       * the filename isn't found then the name will be looked for in the
       * Antiprism data directory colour map resources.
       * \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *map_name, char *errmsg=0);

      ///Initialise from a formatted line
      /**The colour map line is converted to Antiprism format in
       * the following way: ':' is converted to a newline, and '_'
       * to a space. If the lines are bare colours then "idx=" is
       * added.
       * \param map_line the map line text.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the mal line was valid, otherwise false
       * and the error is detailed in \a errmsg. */
      bool init_from_line(const char *map_line, char *errmsg=0);

      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      color_map *clone() { return new color_map_map(*this); }

      ///Get a copy of the map with any with the values indexed by order
      /** \return a pointer to the dynamically allocated condensed copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the copy failed. */
      virtual color_map *get_condensed();

      ///Read the named colours into a colour map
      /**This is equivalent to reading in the x11 colour map from
       * the resources. */
      void read_named_colors();

      ///Set the colour map
      /**\param col_map the colour map to use. */
      void set_map(const map<int, col_val> col_map) { cmap = col_map; }
      
      ///Generate a colour map with HSV values
      /** The map is generated by inserting the H, S, V and A values
       * evenly in map and interpolating the values for the other entries.
       * \param map_sz the number of entries in the map
       * \param hue the list of Hue values
       * \param sat the list of Saturation values
       * \param val the list of Value values
       * \param alpha the list of Alpha values
       * \return \c true if a map was set, otherwise \c false. */
      bool generate_map_hsv(int map_sz, const vector<double> &hue,
            const vector<double> &sat, const vector<double> &val,
            const vector<double> &alpha = vector<double>());
      
      ///Generate a colour map with RGB values
      /** The map is generated by inserting the R, G, V and A values
       * evenly in map and interpolating the values for the other entries.
       * \param map_sz the number of entries in the map
       * \param red the list of Red values
       * \param green the list of Blue values
       * \param blue the list of Green values
       * \param alpha the list of Alpha values
       * \return \c true if a map was set, otherwise \c false. */
      bool generate_map_rgb(int map_sz, const vector<double> &red,
            const vector<double> &green, const vector<double> &blue,
            const vector<double> &alpha = vector<double>());

      ///Generate a colour map
      /** The map name is of the form cmNUMBEROFCOLORS_HUELIST_SATLIST_VALLIST
       * with the lists containing floating point values 0.0 - 1.0 separated
       * by commas e.g. cm256_0,1_0.9_0.9 gives a 256 colour spectrum map.
       * The map is generated by inserting the H, S and V values
       * evenly in map and interpolating the values for the other entries.
       * \param map_name the name of the map
       * \return \c true if a map was set, otherwise \c false. */
      bool generate_map(const char *map_name);

       ///Get the colour map
      /**\return the colour map */
      const map<int, col_val> &get_map() const { return cmap; }
      
      ///Assign operator
      /**\param col_map the colour map to set equal to. */
      color_map &operator =(const color_map &col_map);

      //Reset
      //void reset() { clear(); shift=0; step=1; wrap=0; }
      
      ///Clear the colour map
      void clear() { cmap.clear(); }
      
      ///The size of the map
      /**\return The size */
      unsigned int size() const { return cmap.size(); }
      
      ///The effective size of the map
      /**The effective size of a map is one greater than the highest
       * index number in the map. It is the size of the smallest
       * map (sequential, starting at 0) that will include all
       * the entries of the map.
       * \return The effective size */
      virtual int effective_size() const;
      
      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx);

      ///Set a colour value for an index number.
      /**\param idx the index.
       * \param col the colour to set.
       * \return The colour. */
      void set_col(int idx, col_val col);

      ///Get a colour map containing the colours in this map indexed sequentially.
      /**\param col_map to return the sequential indexing. */
      void get_sequential_index(color_map &col_map);

};

///A colour map that looks up in other colour maps in order
class color_map_multi : public color_map
{
   private:
      //The colour maps to be tried sequentially.
      vector<color_map *> cmaps;

      //Largest effective size
      int max_eff_map_sz;
      
      //Specified map size
      int map_sz;


      //Set map_sz to the largest effective size
      void set_max_eff_map_sz();

   public:
      // Constructor
      color_map_multi() : color_map(), max_eff_map_sz(0) {}

      ///Copy Constructor
      /**\param cmap the multiple colour map to copy from. */
      color_map_multi(const color_map_multi &cmap);
      
      ///Copy Assignment
      /**\param cmap the multiple colour map to copy from. */
      color_map_multi &operator =(const color_map_multi &cmap);
      
      ///Destructor
      ~color_map_multi();

      ///Initialise from a string
      /**The colour map can be in the Antiprism, GIMP or Fractint format. If
       * the filename isn't found then the name will be looked for in the
       * Antiprism data directory colour map resources.
       * \param map_name the map name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool init(const char *map_name, char *errmsg=0);
      
      
      ///Get a copy of the map
      /** \return a pointer to the dynamically allocated copy,
       * which must be freed by the caller with \c delete, 0 indicates
       * that the clone failed. */
      virtual color_map *clone() { return new color_map_multi(*this); }

      ///Set the map size
      /**\param sz the number of entries in the map. */
      void set_map_sz(int sz) { map_sz = sz; }

      ///The effective size of the map
      /**The effective size of a map is one greater than the highest
       * index number in the map. It is the size of the smallest
       * map (sequential, starting at 0) that will include all
       * the entries of the map.
       * \return The effective size */
      virtual int effective_size() const
         { return (map_sz>=0) ? map_sz : max_eff_map_sz; }
      
      ///Add a colour map.
      /**\param col_map the colour map.
       * \param pos the position to add it, or at the end if pos is
       * greater then the current size */
      void add_cmap(color_map *col_map, unsigned int pos=INT_MAX);

      ///Delete a colour map.
      /**\param pos the position of the colour map to delete, or delete
       * the last colour map if \c pos is greater than or equal to the
       * current size */
      void del_cmap(unsigned int pos=INT_MAX);
      
      ///Get the colour value for an index number.
      /**\param idx the index.
       * \return The colour. */
      virtual col_val get_col(int idx);

};

///Create a colour map from its name
/**The map may be read from a file or generated.
 * \param map_name the map name
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return a pointer to the dynamically allocated map which must
 * be freed by the caller with \c delete, if is returned for an
 * invalid map name and the error is detailed in \a errmsg. */
color_map* init_color_map(const char *map_name, char *errmsg=0);




// inline functions

inline int color_map::get_effective_index(int idx) const
{
   int eff_idx = shift + idx*step;
   if(wrap>0)
      eff_idx %= wrap;
   return eff_idx;
}


inline void color_map::cycle_colors()
{
   step++;
   if(wrap)
      step %= wrap;
}


#endif // COL_MAP_H


