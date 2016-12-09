/*
   Copyright (c) 2008-2016, Adrian Rossiter

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

/**\file displaypoly.h
   \brief support for displaying a model as a polyhedron
*/

#ifndef DISPLAYPOLY_H
#define DISPLAYPOLY_H

#include "programopts.h"
#include "scene.h"

namespace anti {

class DisplayPoly : public virtual GeometryDisplay {
private:
  Coloring clrngs[3];
  Timer cmap_tmrs[3];

  bool triangulate;
  unsigned int winding_rule;
  int face_alpha;
  bool use_lines;                    // vrml
  std::vector<std::string> includes; // pov

protected:
  Geometry disp_geom;

  void vrml_trans_begin(FILE *ofile, const Scene &scene);
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

public:
  DisplayPoly();

  Coloring *get_clrngs() { return clrngs; }
  Coloring clrng(int type) { return clrngs[type]; }
  Color def_col(int type);

  void set_triangulate(bool tri);
  bool get_triangulate() { return triangulate; }
  bool set_winding_rule(unsigned int winding);
  unsigned int get_winding_rule() { return winding_rule; }
  void set_face_alpha(int alpha);

  void set_use_lines(bool lines) { use_lines = lines; }
  bool get_use_lines() { return use_lines; }
  void set_includes(std::vector<std::string> incs) { includes = incs; }
  std::vector<std::string> &get_includes() { return includes; }
  const std::vector<std::string> &get_includes() const { return includes; }

  Geometry &get_disp_geom() { return disp_geom; }
  GeometryDisplay *clone() const { return new DisplayPoly(*this); }
  void geom_changed();
  int animate();

  void vrml_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  void pov_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  // void gl_geom(const Scene &scen);
};

class DisplayNumLabels : public virtual GeometryDisplayLabel {
private:
  void vrml_protos(FILE *ofile, const Scene &scen);
  void vrml_verts(FILE *ofile);
  void vrml_edges(FILE *ofile);
  void vrml_faces(FILE *ofile);

  void gl_faces(const Scene &scen);
  void gl_verts(const Scene &scen);
  void gl_edges(const Scene &scen);

public:
  DisplayNumLabels();

  GeometryDisplay *clone() const { return new DisplayNumLabels(*this); };
  void geom_changed() {}

  void vrml_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  void pov_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  // void gl_geom(const scene &scen);
};

class DisplaySymmetry : public virtual DisplayPoly {
protected:
  bool show_axes;
  bool show_mirrors;
  bool show_rotrefls;

  Symmetry sym;
  // Geometry elems;
  // scene_geom sc_geom_elems;

  virtual void disp_changed();

public:
  DisplaySymmetry();

  void set_show_rotrefls(bool show)
  {
    show_rotrefls = show;
    disp_changed();
  }
  bool get_show_rotrefls() { return show_rotrefls; }
  void set_show_mirrors(bool show)
  {
    show_mirrors = show;
    disp_changed();
  }
  bool get_show_mirrors() { return show_mirrors; }
  void set_show_axes(bool show)
  {
    show_axes = show;
    disp_changed();
  }
  bool get_show_axes() { return show_axes; }

  GeometryDisplay *clone() const { return new DisplaySymmetry(*this); }
  void geom_changed();

  void vrml_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  void pov_geom(FILE *ofile, const Scene &scen, int sig_dgts = DEF_SIG_DGTS);
  // void gl_geom(const scene &scen);
};

class ViewOpts : public ProgramOpts {
private:
  DisplayPoly *geom_defs;
  DisplayNumLabels *lab_defs;
  DisplaySymmetry *sym_defs;

public:
  Scene scen_defs;
  Camera cam_defs;

  std::vector<std::string> ifiles;

  static const char *help_view_text;
  static const char *help_scene_text;
  static const char *help_prec_text;

public:
  ViewOpts(const char *name);
  ~ViewOpts();

  Status read_disp_option(char opt, char *optarg);
  void set_view_vals(Scene &scen);
  void set_geom_defs(const DisplayPoly &defs);
  void set_num_label_defs(const DisplayNumLabels &defs);
  void set_sym_defs(const DisplaySymmetry &defs);
  DisplayPoly &get_geom_defs() { return *geom_defs; }
};

} // namespace anti

#endif // DISPLAYPOLY_H
