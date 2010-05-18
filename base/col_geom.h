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

/*!\file col_geom.h
   \brief A class to represent colours in geometries
*/


#ifndef COL_GEOM_H
#define COL_GEOM_H

#include <map>
#include "col_val.h"

using std::map;

///Holder for element colours
class col_geom
{
   private:
      // index of each element type in elem_cols
      enum { v_elems, e_elems, f_elems };

      // Colours for vertices, edges and faces
      map<int, col_val> elem_cols[3];
      
      ///Map the colours to different index numbers.
      /**Used to maintain colors when index numbers are changed. This
       * can happen after deletions.
       * \param cols current element-index-to-colour map.
       * \param chg_map a map of old index numbers to new index numbers.
       *           if the new index number is \c -1 then the element index
       *           has been deleted so the colour is deleted. */
      static void remap_cols(map<int, col_val> &cols,
            const map<int, int> &chg_map);


   protected:
      ///Append a geometry colour holder.
      /**\param geom geometry colour holder to append.
       * \param v_size number of vertices in geometry associated with geom.
       * \param e_size number of edges in geometry associated with geom.
       * \param f_size number of faces in geometry associated with geomi. */
      void append(const col_geom &geom, int v_size,int e_size,int f_size);
     
      ///Map the vertex colours to different index numbers.
      /**Used to maintain colors when index numbers are changed. This
       * can happen after deletions.
       * \param chg_map a map of old index numbers to new index numbers.
       *           if the new index number is \c -1 then the element index
       *           has been deleted so the colour is deleted. */
      void remap_vert_cols(const map<int, int> &chg_map);

      ///Map the edge colours to different index numbers.
      /**Used to maintain colors when index numbers are changed. This
       * can happen after deletions.
       * \param chg_map a map of old index numbers to new index numbers.
       *           if the new index number is \c -1 then the element index
       *           has been deleted so the colour is deleted. */
      void remap_edge_cols(const map<int, int> &chg_map);

      ///Map the face colours to different index numbers.
      /**Used to maintain colors when index numbers are changed. This
       * can happen after deletions.
       * \param chg_map a map of old index numbers to new index numbers.
       *           if the new index number is \c -1 then the element index
       *           has been deleted so the colour is deleted. */
      void remap_face_cols(const map<int, int> &chg_map);

 
   public:
      ///Set a vertex colour.
      /**\param idx the vertex index number.
       * \param col the colour to set. */
      void set_v_col(int idx, col_val col);
      
      ///Get a vertex colour.
      /**\param idx the vertex index number.
       * \return The colour. */
      col_val get_v_col(int idx) const;

      ///Set an edge colour.
      /**\param idx the edge index number.
       * \param col the colour to set. */
      void set_e_col(int idx, col_val col);
      
      ///Get an edge colour.
      /**\param idx the edge index number.
       * \return The colour. */
      col_val get_e_col(int idx) const;

      ///Set a face colour.
      /**\param idx the face index number.
       * \param col the colour to set. */
      void set_f_col(int idx, col_val col);
      
      ///Get a face colour.
      /**\param idx the face index number.
       * \return The colour. */
      col_val get_f_col(int idx) const;

      ///Clear all the vertex colours.
      void clear_v_cols();

      ///Clear all the face colours.
      void clear_e_cols();
     
      ///Clear all the edge colours.
      void clear_f_cols();
      
      ///Clear all the colours.
      void clear_cols();

       ///Read access to vertex colours
      /**\return The colour map */
      const map<int, col_val> &vert_cols() const;
      
      ///Write access to vertex colours
      /**\return The colour map */
      map<int, col_val> &raw_vert_cols();

      ///Read access to edge colours
      /**\return The colour map */
      const map<int, col_val> &edge_cols() const;
      
      ///Write access to edge colours
      /**\return The colour map */
      map<int, col_val> &raw_edge_cols();

      ///Read access to face colours
      /**\return The colour map */
      const map<int, col_val> &face_cols() const;

      ///Write access to face colours
      /**\return The colour map */
      map<int, col_val> &raw_face_cols();

      ///Set a colour in an element-index-to-colour map
      /**\param elem element-index-to-colour map
       * \param idx element index
       * \param col colour to set */
      static void set_col(map<int, col_val> &elem, int idx, col_val col);

      ///Get a colour from an element-index-to-colour map
      /**\param elem element-index-to-colour map
       * \param idx element index to get the colour for
       * \return col The colour. The colour will be in the unset state
       * if the element did not have a colour. */
      static col_val get_col(const map<int, col_val> &elem, int idx);
};




// Implementation of inline functions

inline void col_geom::remap_vert_cols(const map<int, int> &chg_map)
{
   remap_cols(elem_cols[v_elems], chg_map);
}

inline void col_geom::remap_edge_cols(const map<int, int> &chg_map)
{
   remap_cols(elem_cols[e_elems], chg_map);
}

inline void col_geom::remap_face_cols(const map<int, int> &chg_map)
{
   remap_cols(elem_cols[f_elems], chg_map);
}



inline void col_geom::set_v_col(int idx, col_val col)
{
   set_col(elem_cols[v_elems], idx, col);
}

inline col_val col_geom::get_v_col(int idx) const
{
   return get_col(elem_cols[v_elems], idx);
}

inline void col_geom::set_e_col(int idx, col_val col)
{
   set_col(elem_cols[e_elems], idx, col);
}

inline col_val col_geom::get_e_col(int idx) const
{
   return get_col(elem_cols[e_elems], idx);
}

inline void col_geom::set_f_col(int idx, col_val col)
{
   set_col(elem_cols[f_elems], idx, col);
}

inline col_val col_geom::get_f_col(int idx) const
{ 
   return get_col(elem_cols[f_elems], idx);
}

inline void col_geom::clear_v_cols()
{
   elem_cols[v_elems].clear();
}

inline void col_geom::clear_e_cols()
{
   elem_cols[e_elems].clear();
}

inline void col_geom::clear_f_cols()
{
   elem_cols[f_elems].clear();
}

inline void col_geom::clear_cols()
{
   clear_v_cols();
   clear_e_cols();
   clear_f_cols();
}



inline const map<int, col_val> &col_geom::vert_cols() const
{
   return elem_cols[v_elems];
}

inline map<int, col_val> &col_geom::raw_vert_cols()
{
   return elem_cols[v_elems];
}

inline const map<int, col_val> &col_geom::edge_cols() const
{
   return elem_cols[e_elems];
}

inline map<int, col_val> &col_geom::raw_edge_cols()
{
   return elem_cols[e_elems];
}

inline const map<int, col_val> &col_geom::face_cols() const
{
   return elem_cols[f_elems];
}

inline map<int, col_val> &col_geom::raw_face_cols()
{
   return elem_cols[f_elems];
}


inline void col_geom::set_col(map<int, col_val> &elem, int idx, col_val col)
{
   if(col.is_set())
      elem[idx] = col;
   else
      elem.erase(idx);
}  

inline col_val col_geom::get_col(const map<int, col_val> &elem, int idx)
{
   map<int, col_val>::const_iterator mi = elem.find(idx);
   if(mi != elem.end())
      return mi->second;
   else
      return col_val();
}


#endif // COL_GEOM_H

