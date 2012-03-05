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

/*!\file base/geom.h
 * \brief Classes to represent a geometry
*/


#ifndef GEOM_H
#define GEOM_H

#include <vector>
#include <map>
#include <string>

#include "mat3d.h"
#include "vec_utils.h"
#include "col_geom.h"

using std::vector;
using std::map;
using std::string;

class geom_info;

/// Geometry Interface
class geom_if
{
   public:
      /// Destructor
      virtual ~geom_if() {}

      virtual vector<vec3d> *get_verts()=0;
      virtual vector<vector<int> > *get_edges()=0;
      virtual vector<vector<int> > *get_faces()=0;

      /// Check whether geometry is set.
      /**\return \c true if the geometry is set, otherwise \c false. */
      operator bool() const { return verts().size()>0; }

      //-------------------------------------------
      // Element Access
      //-------------------------------------------

      /// Read access to the vertices.
      /**\return A reference to vertex coordinates. */
      virtual const vector<vec3d> &verts() const =0;
      
      /// Read/Write access to the vertices.
      /**\return A reference to vertex coordinates. */
      virtual vector<vec3d> &raw_verts() =0;
      
      /// Read access to a vertex.
      /**\param v_idx index number of the vertex.
       * \return A reference to the vertex coordinates. */
      virtual const vec3d &verts(int v_idx) const =0;
      
      /// Read access to the edges.
      /**\return A reference to the edge data. */
      virtual const vector<vector<int> > &edges() const =0;

      /// Read/Write access to the edges.
      /**\return A reference to the edge data. */
      virtual vector<vector<int> > &raw_edges() =0;
      
      /// Read access to an edge.
      /**\param e_idx index number of the edge.
       * \return A reference to the edge data. */
      virtual const  vector<int> &edges(int e_idx) const =0;

      /// Get the index number of a vertex of an edge. 
      /**\param e_idx edge index number.
       * \param v_no the position the vertex appears in the edge, \c 0 or \c 1
       * \return The vertex index number.
       */
      virtual int edges(int e_idx, int v_no) const =0;
      
      /// Get the coordinates of a vertex of an edge. 
      /**\param e_idx edge index number.
       * \param v_no the position the vertex appears in the edge, \c 0 or \c 1
       * \return The vertex coordinates.
       */
      vec3d edge_v(int e_idx, int v_no) const;
      
      /// Read access to the faces.
      /**\return A reference to the face data. */
      virtual const vector<vector<int> > &faces() const =0;
      
      /// Read/Write access to the faces.
      /**\return A reference to the face data. */
      virtual vector<vector<int> > &raw_faces() =0;
      
      /// Read access to a face.
      /**\param f_idx index number of the face.
       * \return A reference to the face data. */
      virtual const vector<int> &faces(int f_idx) const =0;

      /// Get the index number of a vertex of a face. 
      /**\param f_idx face index number.
       * \param v_no the position the vertex appears in the face,
       * \c 0, \c 1, \c 2, ...
       * \return The vertex index number. */
      virtual int faces(int f_idx, int v_no) const =0;
      
      /// Get the index number of a vertex of a face mod the size of the face 
      /**\param f_idx face index number.
       * \param v_no the position the vertex appears in the face,
       * \c 0, \c 1, \c 2, ...
       * \return The vertex index number. */
      virtual int faces_mod(int f_idx, int v_no) const;

      /// Get the coordinates of a vertex of an face. 
      /**\param f_idx face index number.
       * \param v_no the posisition the vertex appears in the face,
       * \c 0, \c 1, \c 2, ...
       * \return The vertex coordinates.
       */
      vec3d face_v(int f_idx, int v_no) const;



      //-------------------------------------------
      // Add and Delete Elements
      //-------------------------------------------

      /// Add a vertex.
      /**\param vert vertex to add.
       * \return index number of newly added vertex. */
      virtual int add_vert(vec3d vert);
      
      /// Add several vertices.
      /**\param vrts vertices to add.
       * \return index number of last added vertex. */
      virtual int add_verts(const vector<vec3d> &vrts);

      /// Add a face.
      /**\param face face to add.
       * \return index number of newly added face. */
      virtual int add_face(const vector<int> &face);

      /// Add a face made from integer arguments.
      /**The final argument must be a dummy value of \c -1.
       * \param v1 the first index
       * \param ... further indexes in the face, and a final -1 to terminate.
       * \return index number of newly added face. */
      virtual int add_face(int v1, ...);

      /// Add several faces.
      /**\param fces faces to add.
       * \return index number of last added face. */
      virtual int add_faces(const vector<vector<int> > &fces);

      /// Add an edge.
      /**The edge is added only if it is not already present in the edge list.
       * \param edge edge to add.
       * \return index number of edge in edge list.*/
      virtual int add_edge(const vector<int> edge);

      ///Add an edge from vertex index numbers
      /**\param v_idx1 index number of first vertex.
       * \param v_idx2 index number of second vertex.
       * \return index number of edge in edge list.*/
      virtual int add_edge(int v_idx1, int v_idx2);


      /// Add several edges.
      /**Each edge is added only if it is not already present in the edge list.
       * \param edgs edges to add.
       * \return index number of last added edge in edge list. */
      virtual int add_edges(const vector<vector<int> > &edgs);

      /// Add an edge without checking it is valid.
      /**This function is much faster than \c add_edge().
       * The edge should have the lowest vertex index first. The edge
       * should not already be included in the edge list.
       * \param edge edge to add.
       * \return index number of newly added edge. */
      virtual int add_edge_raw(const vector<int> &edge);

      /// Add several edges without checking they are valid.
      /**This function is much faster than \c add_edges().
       * The edges should have the lowest vertex index first. The edges
       * should not already be included in the edge list.
       * \param edgs edges to add.
       * \return index number of last added edge. */
      virtual int add_edges_raw(const vector<vector<int> > &edgs);
      
      /// Delete a vertex
      /**When deleting several vertices is is more efficient to use
       * \c delete_verts().
       * \param v_idx vertex index number to delete
       * \param vert_map a map of old index numbers to new index numbers,
       * deleted vertices map to index \c -1.  */
      virtual void delete_vert(int v_idx, map<int, int> *vert_map=0);
      
      /// Delete several vertices
      /**\param v_idxs vertex index numbers to delete
       * \param vert_map a map of old index numbers to new index numbers,
       * deleted vertices map to index \c -1.  */
      virtual void delete_verts(const vector<int> &v_idxs,
            map<int, int> *vert_map=0);
      
      /// Delete several edges
      /**\param e_idxs edge index numbers to delete
       * \param edge_map a map of old index numbers to new index numbers,
       * deleted edges map to index \c -1.  */
      virtual void delete_edges(const vector<int> &e_idxs, 
            map<int, int> *edge_map=0);
      
      /// Delete several faces
      /**\param f_idxs face index numbers to delete
       * \param face_map a map of old index numbers to new index numbers,
       * deleted faces map to index \c -1.  */
      virtual void delete_faces(const vector<int> &f_idxs,
           map<int, int> *face_map=0);
      
      /// Delete all the vertices
      virtual void clear_verts();
      
      /// Delete all the edges
      virtual void clear_edges();
      
      /// Delete all the faces
      virtual void clear_faces();
      
      /// Delete all the vertices, edges and faces
      virtual void clear_all();
      
      /// Append a geometry
      /**Include the elements of a geometry after the current set of elements*/
      virtual void append(const geom_if &geom);
      
      /// Merge some of the vertices
      /** \param vmap a map from each vertex index number to the index
       * number of the vertex it will be replaced with */
      void verts_merge(map<int, int> &vmap);
      
      /// Add missing implicit edges
      /** Add implicit edges (edges of faces) to the edge list if they are
       * not already included */
      virtual void add_missing_impl_edges();
      
      //-------------------------------------------
      // Transformations
      //-------------------------------------------

      ///Add the faces of the convex hull to the geometry.
      /**\param qh_args additional arguments to pass to qhull (unsupported,
       * may not work, check output.)
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return order or dimension 3, 2, 1 or 0. or -1 if qhull fails
       * and the error is detailed in \a errmsg. */
      int add_hull(string qh_args="", char *errmsg=0);

      ///Set the geometry to its convex hull.
      /**If the convex hull could not be calculated the the geometry
       * will be left empty.
       * \param qh_args additional arguments to pass to qhull (unsupported,
       * may not work, check output.)
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return order or dimension 3, 2, 1 or 0. or -1 if qhull fails
       * and the error is detailed in \a errmsg. */
      int set_hull(string qh_args="", char *errmsg=0);

      ///Get a star of vectors to use for making a zonohedron.
      /**\param type the type of star to make can be
       * <ul>
       * <li>v - centre to vertices are vectors (default)
       * <li>a - all vertex to vertex are vectors
       * <li>i - implicit edges (face sides) are vectors
       * <li>e - explicit edges are vectors
       * </ul>
       * \param centre the centre to use when the type is v.
       * \return The star of vectors. */
      vector<vec3d> get_star(char type='v', vec3d centre=vec3d(0,0,0));

      ///Set the geometry to a zonohedron made from a star
      /**\param star the star of vectors.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the zonohedron could be calculated, otherwise false
       * and the error is detailed in \a errmsg. */
      bool set_zono(const vector<vec3d> &star, char *errmsg=0);

      ///Set planar geodesic division.
      /**A Class I pattern is made with m=0,n=1. A Class II pattern
       * is made with m=1,n=1. A Class III pattern is made with n>=1,m>1.
       * \param base the base polyhedron
       * \param m the first pattern specifier.
       * \param n the second pattern specifier.
       * \return true if the pattern was valid, otherwise false. */
      bool set_geodesic_planar(const geom_if &base, int m, int n=0);

      ///Set spherical geodesic division.
      /**A Class I pattern is made with m=0,n=1. A Class II pattern
       * is made with m=1,n=1. A Class III pattern is made with n>=1,m>1.
       * \param base the base polyhedron
       * \param m the first pattern specifier.
       * \param n the second pattern specifier.
       * \param cent the centre of projection.
       * \return true if the pattern was valid, otherwise false. */
      bool set_geodesic_sphere(const geom_if &base, int m, int n=0,
            vec3d cent=vec3d(0,0,0));

      ///Project the vertices onto a sphere
      /**\param centre the centre of the sphere.
       * \param radius the radius of the sphere. */
      void sphere_projection(vec3d centre=vec3d(0,0,0), double radius=1.0);

      ///Orient the geometry (if possible.)
      /**\param parts used to return the index numbers of the faces
       * in each of the disconnected parts (if \a parts is not \c 0 .)
       * \return the number of disconnected parts. */
      int orient(vector<vector<int> > *parts=0);

      ///Reverse the orientation of all the faces.
      void orient_reverse();

      ///Apply a transformation matrix.
      /**\param trans the transformation matrix. */
      void transform(const mat3d &trans);

      ///Align a polyhedron with the standard alignment for its symmetry type.
      void sym_align();

      ///Triangulate (tesselate) the faces
      /**Divide the faces into triangles, adding new vertices as necessary,
       * and coloured edges if a colour is set for new elements.
       * \param col the colour for any new edges or vertices that were added.
       * The default leaves them with default colours. If it is set to
       * col_val::invisible then the new elements are given a colour that
       * indicates that they should not be displayed.
       * \param winding selects resulting faces acording to winding number
       *    TESS_WINDING_ODD
       *    TESS_WINDING_NONZERO (default)
       *    TESS_WINDING_POSITIVE
       *    TESS_WINDING_NEGATIVE
       *    TESS_WINDING_ABS_GEQ_TWO
       * \param fmap a vector to return the face mapping. Each old face
       * index maps to the first index of faces it was converted to. A
       * final index holds the total number of new faces */
      void triangulate(col_val col=col_val(),
            unsigned int winding=TESS_WINDING_NONZERO, vector<int> *fmap=0);

      
      //-------------------------------------------
      // Geometric Utilities
      //-------------------------------------------

      /// Get the centroid of all the vertices.
      /**\return The coordinates of the centroid. */
      vec3d centroid() const;

      /// Get all the face centroids.
      /**\param ctds vector used to return the coordinates of the centroids. */
      void face_cents(vector<vec3d> &ctds) const;

      /// Get all the face normals.
      /**\param norms vector used to return the normals.
       * \param allow_zero if \c true then the magnitude of the normal
       * is the area of the face, if \c false then for faces with near-zero
       * area a normal with a more accurate direction is calculated. */
      void face_norms(vector<vec3d> &norms, bool allow_zero=false) const;

      /// Get a face centroid, for a face index.
      /**\param f_idx face index number.
       * \return The coordinates of the centroid. */
      vec3d face_cent(int f_idx) const;
      
      /// Get a face centroid, for a face.
      /**\param face contains the vertex index numbers.
       * \return The coordinates of the centroid. */
      vec3d face_cent(const vector<int> &face) const;

      /// Get a face normal, for a face index.
      /**\param f_idx the face index number.
       * \param allow_zero if \c true then the magnitude of the normal
       * is the area of the face, if \c false then for faces with near-zero
       * area a normal with a more accurate direction is calculated.
       * \return The normal */
      vec3d face_norm(int f_idx, bool allow_zero=false) const;
      
      /// Get a face normal, for a face.
      /**\param face contains the vertex index numbers in the face.
       * \param allow_zero if \c true then the magnitude of the normal
       * is the area of the face, if \c false then for faces with near-zero
       * area a normal with a more accurate direction is calculated.
       * \return The normal */
      vec3d face_norm(const vector<int> &face, bool allow_zero=false) const;

      /// Get the nearest point on a face to another point, for a face index.
      /**\param f_idx face index number.
       * \param P the point
       * \return The coordinates of the nearest point to \c P on the
       * face plane. */
      vec3d face_nearpt(int f_idx, vec3d P) const;

      /// Get the nearest point on a face to another point, for a face.
      /**\param face contains the vertex index numbers.
       * \param P the point
       * \return The coordinates of the nearest point to \c P on the
       * face plane. */
      vec3d face_nearpt(const vector<int> &face, vec3d P) const;
      
      /// Get an edge centre (the centroid), for an edge index.
      /**\param e_idx edge index number.
       * \return The coordinates of the centre. */
      vec3d edge_cent(int e_idx) const;

      /// Get an edge centre (the centroid), for an edge.
      /**\param edge contains the two vertex index numbers.
       * \return The coordinates of the centre. */
      vec3d edge_cent(const vector<int> &edge) const;
      
      /// Get the vector between the ends of an edge, for an edge index.
      /**\param e_idx edge index number.
       * \return The vector. */
      vec3d edge_vec(int e_idx) const;
      
      /// Get the vector between the ends of an edge, for an edge.
      /**\param edge contains the two vertex index numbers.
       * \return The vector. */
      vec3d edge_vec(const vector<int> &edge) const;
      
      /// Get the vector between two vertices by vertex index.
      /**\param v_idx0 start vertex index number.
       * \param v_idx1 end vertex index number.
       * \return The vector. */
      vec3d edge_vec(int v_idx0, int v_idx1) const;

      /// Get the nearest point on an edge to another point, for a edge index.
      /**\param e_idx edge index number.
       * \param P the point
       * \return The coordinates of the nearest point to \c P on the
       * line of the edge. */
      vec3d edge_nearpt(int e_idx, vec3d P) const;

      /// Get the nearest point on an edge to another point, for an edge.
      /**\param edge contains the two vertex index numbers.
       * \param P the point
       * \return The coordinates of the nearest point to \c P on the
       * line of the edge. */
      vec3d edge_nearpt(const vector<int> &edge, vec3d P) const;
      
      /// Get the length of an edge, for an edge index.
      /**\param e_idx edge index number.
       * \return The length of the edge. */
      double edge_len(int e_idx) const;
      
      /// Get the length of an edge, for an edge.
      /**\param edge contains the two vertex index numbers.
       * \return The length of the edge. */
      double edge_len(const vector<int> &edge) const;


      //-------------------------------------------
      //  Other Utilities
      //-------------------------------------------
     
      ///Read geometry from a file
      /**The file is first read as a normal OFF file, if that fails it will be
       * read as a Qhull formatted OFF file, and if that fails the file will be
       * read for any coordinates (lines that contains three numbers separated
       * by commas and/or spaces will be taken as a set of coordinates.)
       * \param file_name the file name ("" or "-" for standard input).
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool read(string file_name="", char *errmsg=0);

      ///Read geometry from a file stream
      /**The file is first read as a normal OFF file, if that fails it will be
       * read as a Qhull formatted OFF file, and if that fails the file will be
       * read for any coordinates (lines that contains three numbers separated
       * by commas and/or spaces will be taken as a set of coordinates.)
       * \param file the file stream.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the file could be read, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool read(FILE *file, char *errmsg=0);

      ///Read resource geometry from resource model name
      /**\param res_name the resource model name.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if the resource name was valid, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool read_resource(string res_name="", char *errmsg=0);

      ///Write geometry to a file
      /**\param file_name the file name ("" for standard output.)
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \param sig_dgts the number of significant digits to write,
       * or if negative then the number of digits after the decimal point.
       * \return true if the file could be written, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool write(string file_name="", char *errmsg=0,
            int sig_dgts=DEF_SIG_DGTS) const;

      ///Write geometry to a file stream
      /**\param file the file stream.
       * \param sig_dgts the number of significant digits to write,
       * or if negative then the number of digits after the decimal point. */
      virtual void write(FILE *file, int sig_dgts=DEF_SIG_DGTS) const;

      ///Write coordinates to a file
      /**\param file_name the file name ("" for standard output.)
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \param sep a string to use as the seperator between coordinates.
       * \param sig_dgts the number of significant digits to write,
       * or if negative then the number of digits after the decimal point.
       * \return true if the file could be written, otherwise false
       * and the error is detailed in \a errmsg. */
      virtual bool write_crds(string file_name="", char *errmsg=0,
            const char *sep=" ", int sig_dgts=DEF_SIG_DGTS) const;

      ///Write coordinates to a file stream
      /**\param file the file stream.
       * \param sep a string to use as the seperator between coordinates.
       * \param sig_dgts the number of significant digits to write,
       * or if negative then the number of digits after the decimal point. */
      virtual void write_crds(FILE *file, const char *sep=" ",
            int sig_dgts=DEF_SIG_DGTS) const;

      ///Check if geomtery is consistently oriented
      /**\return \c true if consistently oriented, otherwise \c false. */
      bool is_oriented() const;

      /// Get geom_info object
      /**\return geom_info object associated with this geometry. */
      geom_info get_info() const;
      
      /// Get implicit edges
      /** Returns the edges of the polygon faces
       * \param edgs the edges are returned here */
      void get_impl_edges(vector<vector<int> > &edgs) const;
      void get_edge_face_pairs(map<vector<int>, vector<int> > &edge2facepr,
            bool oriented=true) const;
      
};

///Make an edge from vertex index numbers
/**\param v_idx1 index number of first vertex.
 * \param v_idx2 index number of second vertex.
 * \return edge containg index numbers in numerical order. */
vector<int> make_edge(int v_idx1, int v_idx2);


/// Geometry that holds the elements by value.
class geom_v: public geom_if
{
   private:
      vector<vec3d> vert_elems;
      vector<vector<int> > face_elems;
      vector<vector<int> > edge_elems;

   public:
      /// Constructor
      geom_v() {}

      /// Constructor
      /** Initialise from another geometry that implements \c geom_if */
      geom_v(const geom_if &geom)
        :  vert_elems(geom.verts()),
           face_elems(geom.faces()),
           edge_elems(geom.edges())
           {}

      vector<vec3d> *get_verts() {return &vert_elems; }
      vector<vector<int> > *get_edges() { return &edge_elems; }
      vector<vector<int> > *get_faces() { return &face_elems; }
      
      virtual const vector<vec3d> &verts() const;
      virtual vector<vec3d> &raw_verts();
      virtual const vec3d &verts(int v_idx) const;
      
      virtual const vector<vector<int> > &edges() const;
      virtual vector<vector<int> > &raw_edges();
      virtual const  vector<int> &edges(int e_idx) const;
      virtual int edges(int e_idx, int v_no) const;
      
      virtual const vector<vector<int> > &faces() const;
      virtual vector<vector<int> > &raw_faces();
      virtual const vector<int> &faces(int f_idx) const;
      virtual int faces(int f_idx, int v_no) const;


};


/// Geometry with coloured elements.
class col_geom_v: public col_geom, public geom_v {
   public:
      /// Constructor
      col_geom_v() {}

      /// Copy Constructor
      /** Initialise from another geometry that implements \c geom_if */
      col_geom_v(const geom_if &geom);

      /// Copy Assignment 
      /** Initialise from another geometry that implements \c geom_if */
      col_geom_v &operator =(const geom_if &geom);

      /// Add a vertex with a colour.
      /**\param vert vertex to add.
       * \param col colour of the vertex.
       * \return index number of newly added vertex. */
      int add_col_vert(vec3d vert, col_val col);

      /// Add a face with a colour.
      /**\param face face to add.
       * \param col colour of the face.
       * \return index number of newly added face. */
      int add_col_face(const vector<int> &face, col_val col);

      /// Add an edge with a colour.
      /**\param edge edge to add.
       * \param col colour of the edge.
       * \return index number of newly added edge. */
      int add_col_edge(const vector<int> &edge, col_val col);

      ///Add an edge from vertex index numbers, with a colour.
      /**\param v_idx1 index number of first vertex.
       * \param v_idx2 index number of second vertex.
       * \param col colour of the edge.
       * \return index number of newly added edge.*/
      virtual int add_col_edge(int v_idx1, int v_idx2, col_val col);

      virtual void add_missing_impl_edges();
      
      /// Color convenience function
      virtual void color_vef(col_val vert_col, col_val edge_col, col_val face_col);
       
      virtual void delete_verts(const vector<int> &v_nos,
            map<int, int> *vert_map=0);
      virtual void delete_edges(const vector<int> &e_nos,
            map<int, int> *edge_map=0);
      virtual void delete_faces(const vector<int> &f_nos,
            map<int, int> *face_map=0);
      virtual void append(const geom_if& geom);
      
      virtual void clear_verts();
      virtual void clear_edges();
      virtual void clear_faces();

};






// Implementation of inline functions

inline vec3d geom_if::edge_v(int e_idx, int v_no) const
{
   return verts(edges(e_idx, v_no));
}

inline vec3d geom_if::face_v(int f_idx, int v_no) const
{
   return verts(faces(f_idx, v_no));
}

inline int geom_if::faces_mod(int f_idx, int v_no) const
{
   unsigned int f_sz = faces(f_idx).size();
   v_no = v_no%f_sz;
   if(v_no<0)
      v_no += f_sz;
   return faces(f_idx, v_no);
}



inline vector<vec3d> &geom_v::raw_verts()
{
   return vert_elems;
}

inline const vector<vec3d> &geom_v::verts() const
{
   return vert_elems;
}

inline const vec3d &geom_v::verts(int v_idx) const
{
   return vert_elems[v_idx];
}

inline vector<vector<int> > &geom_v::raw_edges()
{
   return edge_elems;
}

inline const vector<vector<int> > &geom_v::edges() const
{
   return edge_elems;
}

inline const vector<int> &geom_v::edges(int e_idx) const
{
   return edge_elems[e_idx];
}

inline int geom_v::edges(int e_idx, int v_no) const
{
   return edge_elems[e_idx][v_no];
}



inline vector<vector<int> > &geom_v::raw_faces()
{
   return face_elems;
}

inline const vector<vector<int> > &geom_v::faces() const
{
   return face_elems;
}

inline const vector<int> &geom_v::faces(int f_idx) const
{
   return face_elems[f_idx];
}

inline int geom_v::faces(int f_idx, int v_no) const
{
   return face_elems[f_idx][v_no];
}

// col_geom_v

inline int col_geom_v::add_col_vert(vec3d vert, col_val col)
{ 
   int idx = geom_if::add_vert(vert);
   set_v_col(idx, col);
   return idx;
}

inline int col_geom_v::add_col_edge(const vector<int> &edge, col_val col)
{ 
   int idx = geom_if::add_edge(edge);
   set_e_col(idx, col);
   return idx;
}

inline int col_geom_v::add_col_edge(int v_idx1, int v_idx2, col_val col)
{
   return add_col_edge(make_edge(v_idx1, v_idx2), col);
}



inline int col_geom_v::add_col_face(const vector<int> &face, col_val col)
{ 
   int idx = geom_if::add_face(face);
   set_f_col(idx, col);
   return idx;
}


inline void col_geom_v::clear_verts()
{
   geom_if::clear_verts();
   col_geom::clear_v_cols();
}

inline void col_geom_v::clear_edges()
{
   geom_if::clear_edges();
   col_geom::clear_e_cols();
}

inline void col_geom_v::clear_faces()
{
   geom_if::clear_faces();
   col_geom::clear_f_cols();
}



inline vec3d geom_if::centroid() const
{
   return ::centroid(verts());
}

inline void geom_if::face_cents(vector<vec3d> &ctds) const
{
   ctds.resize(faces().size());
   for(unsigned int i=0; i<faces().size(); i++)
      ctds[i] = face_cent(i);
}


inline void geom_if::face_norms(vector<vec3d> &norms, bool allow_zero) const
{
   norms.resize(faces().size());
   for(unsigned int i=0; i<faces().size(); i++)
      norms[i] = face_norm(i, allow_zero).unit();
}


// Faces

inline vec3d geom_if::face_cent(int f_idx) const
{
   return face_cent(faces(f_idx));
}

inline vec3d geom_if::face_cent(const vector<int> &face) const
{
   return ::centroid(verts(), face);
}

inline vec3d geom_if::face_norm(int f_idx, bool allow_zero) const
{
   return face_norm(faces(f_idx), allow_zero);
}

inline vec3d geom_if::face_norm(const vector<int> &face,
      bool allow_zero) const
{
   return ::face_norm(verts(), face, allow_zero);
}

inline vec3d geom_if::face_nearpt(int f_idx, vec3d P) const
{
   return face_nearpt(faces(f_idx), P);
}

inline vec3d geom_if::face_nearpt(const vector<int> &face, vec3d P) const
{
   return nearest_point(P, verts(), face);
}


// Edges

inline vec3d geom_if::edge_cent(int e_idx) const
{
   return edge_cent(edges(e_idx));
}

inline vec3d geom_if::edge_cent(const vector<int> &edge) const
{
   return 0.5 * ( verts(edge[1]) + verts(edge[0]) );
}

inline vec3d geom_if::edge_vec(int e_idx) const
{
   return edge_vec(edges(e_idx));
}

inline vec3d geom_if::edge_vec(const vector<int> &edge) const
{
   return verts(edge[1]) - verts(edge[0]);
}

inline vec3d geom_if::edge_vec(int v_idx0, int v_idx1) const
{
   return verts(v_idx1) - verts(v_idx0);
}

inline vec3d geom_if::edge_nearpt(int e_idx, vec3d P) const
{
   return edge_nearpt(edges(e_idx), P);
}

inline vec3d geom_if::edge_nearpt(const vector<int> &edge, vec3d P) const
{
   return nearest_point(P, verts(), edge);
}

inline double geom_if::edge_len(int e_idx) const
{
   return edge_vec(e_idx).mag();
}

inline double geom_if::edge_len(const vector<int> &edge) const
{
   return edge_vec(edge).mag();
}

inline void geom_if::transform(const mat3d &trans)
{
   ::transform(raw_verts(), trans);
}

inline col_geom_v::col_geom_v(const geom_if &geom)
{ 
   append(geom);
}


#endif // GEOM_H

