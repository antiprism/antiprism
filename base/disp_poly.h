/*
   Copyright (c) 2008, Adrian Rossiter

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

/**\file disp_poly.h
   \brief support for displaying a model as a polyhedron
*/


#ifndef DISP_POLY_H
#define DISP_POLY_H

#include "scene.h"

class disp_poly : public geom_disp
{
   private:
      coloring clrngs[3];
      timer cmap_tmrs[3];

      bool triangulate;
      int  face_alpha;
      bool use_lines;              // vrml
      vector<string> includes;     // pov

   protected:
      col_geom_v disp_geom;

      void vrml_trans_begin(FILE *ofile, const scene &scene);
      void vrml_trans_end(FILE *ofile);
      void vrml_protos(FILE *ofile);
      void vrml_coords(FILE *ofile, int sig_digits);
      void vrml_verts_l(FILE *ofile);
      void vrml_verts(FILE *ofile, int sig_digits);
      void vrml_edges_l(FILE *ofile);
      void vrml_edges(FILE *ofile);
      void vrml_faces(FILE *ofile);

      void pov_default_vals(FILE *ofile);
      void pov_disp_macros(FILE *ofile);
      void pov_vert_arrays(FILE *ofile, int sig_digits);
      void pov_edge_arrays(FILE *ofile);
      void pov_face_arrays(FILE *ofile);
      void pov_elements(FILE *ofile, int sig_digits);
      void pov_col_maps(FILE *ofile);
      void pov_include_files(FILE *ofile);
      void pov_object(FILE *ofile);

      void gl_verts(const scene &scen);
      void gl_edges(const scene &scen);
      void gl_faces(const scene &scen);


   public:
      disp_poly();
      
      coloring *get_clrngs() { return clrngs; }
      coloring &get_v_clrng() { return clrngs[0]; }
      coloring &get_e_clrng() { return clrngs[1]; }
      coloring &get_f_clrng() { return clrngs[2]; }
      
      col_val get_def_v_col();
      col_val get_def_e_col();
      col_val get_def_f_col();
      
      void set_triangulate(bool tri);
      bool get_triangulate() { return triangulate; }
      void set_face_alpha(int alpha);

      void set_use_lines(bool lines) { use_lines = lines; }
      bool get_use_lines() { return use_lines; }
      void set_includes(vector<string> incs) { includes = incs; }
      vector<string> &get_includes() { return includes; }
      const vector<string> &get_includes() const { return includes; }
 
      geom_disp *clone() const { return new disp_poly(*this); }
      void geom_changed();
      void animate();
      
      void vrml_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void pov_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void gl_geom(const scene &scen);
};





class disp_num_labels : public geom_disp_label
{
   private:

      void vrml_protos(FILE *ofile, const scene &scene);
      void vrml_verts(FILE *ofile);
      void vrml_edges(FILE *ofile);
      void vrml_faces(FILE *ofile);
      
      void gl_faces(const scene &scen);
      void gl_verts(const scene &scen);
      void gl_edges(const scene &scen);

   public:
      disp_num_labels();

      geom_disp *clone() const { return new disp_num_labels(*this); };
      void geom_changed() {}

      void vrml_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void pov_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void gl_geom(const scene &scen);
 };


class disp_sym: public disp_poly
{
   protected:
      bool show_axes;
      bool show_mirrors;
      bool show_rotrefls;

      sch_sym sym;
      //col_geom_v elems;
      //scene_geom sc_geom_elems;
      
      virtual void disp_changed();

   public:
      disp_sym();
      
      void set_show_rotrefls(bool show) { show_rotrefls=show; disp_changed();}
      bool get_show_rotrefls() { return show_rotrefls; }
      void set_show_mirrors(bool show) { show_mirrors=show; disp_changed();}
      bool get_show_mirrors() { return show_mirrors; }
      void set_show_axes(bool show) { show_axes=show; disp_changed();}
      bool get_show_axes() { return show_axes; }
      
      geom_disp *clone() const { return new disp_sym(*this); }
      void geom_changed();

      void vrml_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void pov_geom(FILE *ofile, const scene &scen, int sig_dgts=DEF_SIG_DGTS);
      void gl_geom(const scene &scen);

};

class view_opts: public prog_opts {
   public:
      scene scen_defs;
      disp_poly geom_defs;
      disp_num_labels lab_defs;
      disp_sym sym_defs;
      camera cam_defs;
      
      vector<string> ifiles;

      static const char *help_view_text;
      static const char *help_scene_text;
      static const char *help_prec_text;

   public:
      view_opts(const char *name): prog_opts(name) {}
      bool read_disp_option(char opt, char *optarg, char *errmsg);
      void set_view_vals(scene &scen);
};


#endif // DISP_POLY_H


