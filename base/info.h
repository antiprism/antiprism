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
   Name: info.h
   Description: information from OFF files
   Project: Antiprism - http://www.antiprism.com
*/


#ifndef INFO_H
#define INFO_H

#include "geom.h"
#include "transforms.h"
using std::pair;


int cmp_angles(const double &a, const double &b, double min_diff=1e-6);
int cmp_face_angles(const vector<double> &f1, const vector<double> &f2, double min_diff=1e-6);

class ang_less {
   public:
      bool operator() (const double &a, const double &b) {
         return (cmp_angles(a, b) < 0);
      }
};

class ang_vect_less {
   public:
      bool operator() (const vector<double> &f1, const vector<double> &f2) {
         if(f1.size() != f2.size())
            return (f1.size() < f2.size());
         return (cmp_face_angles(f1, f2) < 0);
      }
};


class elem_lims
{
   public:
      enum { IDX_MIN=0, IDX_MIN2, IDX_MAX, IDX_MAX2, IDX_MID, IDX_MID2 };
      
      int idx[6];
      double max;
      double min;
      double mid;
      double sum;
      
      elem_lims() { init(); }
      void init() { idx[0]=-1; max=-1e100; min=1e100; mid=1e100; sum=0; }
      bool is_set() const { return idx[0] != -1; }
};
   
class elem_distances : public elem_lims
{
   protected:
      const geom_if &geom;
      vec3d center;
   public:

      elem_distances(const geom_if &geo): geom(geo) { set_center();}
      virtual ~elem_distances() {}
      void set_center(vec3d cent=vec3d(0,0,0)) {center = cent; idx[0] = -1; }
      virtual void set_values() = 0;
};

class v_distances: public elem_distances
{ 
   public:
      v_distances(const geom_if &geom): elem_distances(geom) { idx[0] = -1; }
      void set_values();
};

class e_distances: public elem_distances
{ 
   public:
      e_distances(const geom_if &geom): elem_distances(geom) { idx[0] = -1; }
      void set_values();
};

class ie_distances: public elem_distances
{
   public:
      ie_distances(const geom_if &geom): elem_distances(geom) { idx[0] = -1; }
      void set_values();
};

class f_distances: public elem_distances
{ 
   public:
      f_distances(const geom_if &geom): elem_distances(geom) { idx[0] = -1; }
      void set_values();
};



class geom_info
{
   private:
      vec3d cent;
      bool closed;
      bool polyhedron;
      int oriented;
      int orientable;
      int number_parts;
      elem_lims iedge_len;
      elem_lims edge_len;
      elem_lims so_angles;
      elem_lims dih_angles;
      elem_lims ang;
      int num_angs;
      elem_lims area;
      double vol;
      v_distances v_dsts;
      e_distances e_dsts;
      ie_distances ie_dsts;
      f_distances f_dsts;

      vector<vector<int> > impl_edges;
      map<vector<int>, vector<int> > efpairs;
      map<vector<double>, int, ang_vect_less> face_angles;
      map<vector<double>, int, ang_vect_less> vert_dihed;
      map<double, int, ang_less> dihedral_angles;
      map<double, int, ang_less> e_lengths;
      map<double, int, ang_less> ie_lengths;
      map<double, int, ang_less> plane_angles;
      map<double, int, ang_less> sol_angles;
      vector<double> vertex_angles;
      map<pair<int, int>, double> vertex_plane_angs;
      vector<double> edge_dihedrals;
      vector<double> f_areas;
      vector<double> f_max_nonplanars;
      vector<vector<int> > vert_cons;
      vector<int> free_verts;
      bool free_verts_found;
      geom_v dual;
      sch_sym sym;
   
      void find_impl_edges() { geom.get_impl_edges(impl_edges); }
      void find_edge_face_pairs();
      void find_face_angles();
      void find_dihedral_angles();
      void find_vert_cons();
      void find_free_verts();
      void find_solid_angles();
      void find_e_lengths(map<double, int, ang_less> &e_lens,
         const vector<vector<int> > &edges, elem_lims &lens);
      void find_e_lengths();
      void find_f_areas();
      void find_f_max_nonplanars();
      void find_oriented();
      void find_symmetry();

   protected:   
      const geom_if &geom;

   public:
      geom_info(const geom_if &geo) :
         cent(vec3d(0,0,0)), oriented(-1), orientable(-1),
         v_dsts(geo), e_dsts(geo), ie_dsts(geo), f_dsts(geo),
         free_verts_found(false), geom(geo)
         {}
      
      void reset();
      void set_center(vec3d center) { cent = center; v_dsts.set_center(cent);
         e_dsts.set_center(cent); ie_dsts.set_center(cent);
         f_dsts.set_center(cent); }

      // elements
      double face_area(int f_no);
      double face_vol(int f_no);
      void face_angles_lengths(int f_no, vector<double> &angs,
            vector<double> *lens=0);
      
      bool is_closed();
      bool is_oriented();
      bool is_orientable();
      bool is_polyhedron();
      
      int num_verts() { return geom.verts().size(); }
      int num_edges() { return geom.edges().size(); }
      int num_iedges() { return get_impl_edges().size(); }
      int num_faces() { return geom.faces().size(); }
      int num_parts() { is_orientable(); return number_parts; }

      elem_lims face_areas() { if(f_areas.size()==0) find_f_areas();
                           return area; }
      double volume() { if(f_areas.size()==0) find_f_areas(); return vol;}
      double vol2_by_area3() { return pow(volume(),2)/pow(face_areas().sum,3); }

      elem_lims edge_lengths() { get_e_lengths(); return edge_len; }
      elem_lims iedge_lengths() { get_ie_lengths(); return iedge_len; }
      elem_lims dihed_angles()
                  { if(dihedral_angles.size()==0) find_dihedral_angles();
                    return dih_angles; }
      elem_lims solid_angles() {if(sol_angles.size()==0) find_solid_angles();
                                return so_angles; }
      vec3d center() { return cent; }
      elem_distances &vert_dists()
         { if(!v_dsts.is_set()) v_dsts.set_values(); return v_dsts; }
      elem_distances &edge_dists()
         { if(!e_dsts.is_set()) e_dsts.set_values(); return e_dsts; }
      elem_distances &impl_edge_dists()
         { if(!ie_dsts.is_set()) ie_dsts.set_values(); return ie_dsts; }
      elem_distances &face_dists()
         { if(!f_dsts.is_set()) f_dsts.set_values(); return f_dsts; }
      elem_lims angles()
         { if(!plane_angles.size()) find_face_angles(); return ang; }
      int num_angles()
         { if(!plane_angles.size()) find_face_angles(); return num_angs; }
      double angle_defect() { return num_verts()*2*M_PI - angles().sum; }

      //verts
      const vector<vector<int> > &get_vert_cons()
         { if(!vert_cons.size()) find_vert_cons(); return vert_cons; }
      const vector<double> &get_vertex_angles()
         { if(!vertex_angles.size()) find_solid_angles(); return vertex_angles;}
      map<double, int, ang_less> &get_solid_angles()
         { if(!sol_angles.size()) find_solid_angles(); return sol_angles; }
      map<pair<int, int>, double> &get_vertex_plane_angs()
         { if(!vertex_plane_angs.size()) find_face_angles();
            return vertex_plane_angs; }
      const vector<int> &get_free_verts()
         { if(!free_verts_found) find_free_verts(); return free_verts;}
      
      //edges
      map<vector<int>, vector<int> > &get_edge_face_pairs()
         { if(!efpairs.size()) find_edge_face_pairs(); return efpairs; }
      vector<double> &get_edge_dihedrals()
         { if(!dihedral_angles.size()) find_dihedral_angles();
            return edge_dihedrals; }
      map<double, int, ang_less> &get_dihedral_angles()
         { if(!dihedral_angles.size()) find_dihedral_angles();
            return dihedral_angles; }
      map<double, int, ang_less> &get_e_lengths()
         { if(!e_lengths.size())
              find_e_lengths(e_lengths, geom.edges(), edge_len);
           return e_lengths; }
      
      //implicit edges
      vector<vector<int> > &get_impl_edges()
         { if(!impl_edges.size()) find_impl_edges(); return impl_edges; }
      map<double, int, ang_less> &get_ie_lengths()
         { if(!ie_lengths.size())
              find_e_lengths(ie_lengths, get_impl_edges(), iedge_len);
           return ie_lengths; }
      
      //faces
      map<vector<double>, int, ang_vect_less> &get_face_angles()
         { if(!face_angles.size()) find_face_angles(); return face_angles; }
      vector<double> &get_f_areas()
         { if(!f_areas.size()) find_f_areas(); return f_areas; }
      vector<double> &get_f_max_nonplanars()
         { if(!f_max_nonplanars.size()) find_f_max_nonplanars();
           return f_max_nonplanars; }
      geom_v &get_dual()
         { if(!dual.get_faces()->size()) ::get_dual(geom, dual); return dual; } 

      const geom_if &get_geom() const { return geom; }
  
      //Symmetry
      string get_symmetry_type_name()
         { if(!sym.is_set()) find_symmetry(); return sym.get_symbol(); }
      const set<sch_axis> &get_symmetry_axes()
         { if(!sym.is_set()) find_symmetry(); return sym.get_axes(); }
      mat3d get_symmetry_alignment_to_std()
         { if(!sym.is_set()) find_symmetry(); return sym.get_to_std(); }
};


#endif // INFO_H
