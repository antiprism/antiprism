
%module antiprism

%{
#include "../base/antiprism.h"
%}

%include "antiprism_doc.i"

%include "std_map.i"
%include "std_set.i"
%include "std_string.i"
%include "std_vector.i"
%include "cstring.i"

%cstring_bounded_mutable(char *errmsg, 256);
%cstring_bounded_mutable(char *err_msg, 256);


%include "../base/antiprism.h"
%include "../base/const.h"

%inline {
   col_geom_v read_geom(string file_name, char *err_msg)
      { col_geom_v g; g.read(file_name, err_msg); return g; }
   col_geom_v read_geom(char *err_msg)
      { return read_geom("", err_msg); }
   vec3d read_vec3d(char *str, char *err_msg)
      { vec3d v; v.read(str, err_msg); return v; }
   vec4d read_vec4d(char *str, char *err_msg)
      { vec4d v; v.read(str, err_msg); return v; }
   col_val read_color(char * str, char *err_msg)
      { col_val col; col.read(str, err_msg); return col; }

}



%rename(__iadd__) vec3d::operator+=;
%rename(__isub__) vec3d::operator-=;
%rename(__imul__) vec3d::operator*=;
%rename(__idiv__) vec3d::operator/=;

%ignore vec3d::read(char *);
%ignore vec3d::operator[];
%ignore operator +(vec3d, vec3d);
%ignore operator -(vec3d, vec3d);
%ignore operator-(vec3d);
%ignore operator -(vec3d, vec3d);
%ignore operator *(vec3d, double);
%ignore operator *(double, vec3d);
%ignore operator /(vec3d, double);
%ignore operator /(double, vec3d);

%extend vec3d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      vec3d operator +(vec3d v) { return *$self + v; }
      vec3d operator -(vec3d v) { return *$self - v; }
      vec3d operator -() { return -(*$self); }
      vec3d operator *(double n) { return *$self * n; }
      vec3d operator /(double n) { return *$self / n; }
}

%include "../base/vec3d.h"



%rename(__neg__) vec4d::operator-();
%rename(__iadd__) vec4d::operator+=;
%rename(__isub__) vec4d::operator-=;
%rename(__imul__) vec4d::operator*=;
%rename(__idiv__) vec4d::operator/=;

%ignore vec4d::read(char *);
%ignore vec4d::operator[];
%ignore operator +(vec4d, vec4d);
%ignore operator -(vec4d, vec4d);
%ignore operator-(vec4d);
%ignore operator -(vec4d, vec4d);
%ignore operator *(vec4d, double);
%ignore operator *(double, vec4d);
%ignore operator /(vec4d, double);
%ignore operator /(double, vec4d);

%extend vec4d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      vec4d operator +(vec4d v) { return *$self + v; }
      vec4d operator -(vec4d v) { return *$self - v; }
      vec4d operator -() { return -(*$self); }
      vec4d operator *(double n) { return *$self * n; }
      vec4d operator /(double n) { return *$self / n; }
}

%include "../base/vec4d.h"



%rename(__imul__) mat3d::operator*=;

%ignore mat3d::operator[];
%ignore operator *(const mat3d &, const vec3d &);
%ignore operator *(const vec3d &, const mat3d &);
%ignore operator *(const mat3d &m1, const mat3d &m2);
%ignore operator <(const mat3d &m1, const mat3d &m2);

%extend mat3d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      vec3d operator *(const vec3d &vec) { return *$self * vec; }
      //vec3d operator *(const mat3d &mat) { return mat * (*$self); }
      mat3d operator *(const mat3d &mat) { return *$self * mat; }
      bool operator <(const mat3d &mat) { return *$self < mat; }
}

%include "../base/mat3d.h"



%ignore line_plane_intersect(vec3d Q, vec3d n, vec3d P0, vec3d P1, int *);

%include "../base/vec_utils.h"



%ignore mat4d::operator[];
%ignore operator *(const mat4d &, const vec4d &);
%ignore operator *(const vec4d &, const mat4d &);
%ignore operator *(const mat4d &, const mat4d &);
%ignore operator *(double, const mat4d &);
%ignore operator *(const mat4d &, double);
%ignore operator +(const mat4d &, const mat4d &);
%rename(__imul__) mat4d::operator*=;
%rename(__iadd__) mat4d::operator+=;

%extend mat4d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      vec4d operator *(const vec4d &vec) { return *$self * vec; }
      //vec4d operator *(const mat4d &mat) { return *$self * mat; }
      mat4d operator *(const mat4d &mat) { return *$self * mat; }
      mat4d operator *(double n) { return *$self * n; }
      mat4d operator +(const mat4d &mat) { return *$self + mat; }
}

%include "../base/mat4d.h"



%ignore col_val::operator[];
%ignore col_val::read(char *);
%ignore col_val::from_offvals(vector<char *> &);
%ignore col_val::from_decvals(vector<double> &);
%ignore col_val::from_intvals(vector<int> &);
%ignore col_val::read_decvals(vector<char *> &);
%ignore col_val::read_intvals(vector<char *> &);
%ignore col_val::read_decvals(char *);
%ignore col_val::read_intvals(char *);
%ignore col_val::read_hexvals(char *);
%ignore col_val::read_colorname(char *);

%extend col_val {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
}

%include "../base/col_val.h"



%include "../base/col_geom.h"



%rename(is_set) geom_if::operator bool() const;
%rename(__eq__) col_geom_v::operator=(const geom_if &);

%ignore geom_if::add_hull(string);
%ignore geom_if::set_hull(string);
%ignore geom_if::read(string);
%ignore geom_if::write(string);
%ignore geom_if::write_crds(string);
%ignore geom_if::read(FILE *);
%ignore geom_if::read(FILE *, char *);
%ignore geom_if::write(FILE *);
%ignore geom_if::write(FILE *, char *);
%ignore geom_if::write(FILE *, char *, int);

%extend geom_if {
   //write() { return $self->write(""); }

}

%include "../base/geom.h"



%rename(__eq__) color_map::operator=;

%ignore color_map::read(string);

%include "../base/col_map.h"



%include "../base/coloring.h"



%rename(is_set) t_set::operator bool() const;

%ignore operator *(const t_set &, const t_set &);
%ignore operator *=(t_set &, const t_set &);
%ignore operator +(const t_set &, const mat3d &);
%ignore operator +=(t_set &, const mat3d &);
%ignore sch::read(string);

%extend t_set {
   public:
      t_set operator *(const t_set &s) { return t_set().product(*$self, s); }
      t_set &operator *=(const t_set &s) { return $self->product_with(s); }
      t_set operator +(const mat3d &m) { return t_set(*$self).add(m); }
      t_set &operator +=(const mat3d &m) { return $self->add(m); }
}

t_set &operator *=(t_set &s1, const t_set &s2);

%include "../base/symmetry.h"


%rename(__eq__) scene_geom::operator=(const scene_geom &);

%rename(__eq__) scene::operator=(const scene &);


%include "../base/scene.h"


%include "../base/transforms.h"



%ignore minimum_distance(const geom_if &geom, double);

%include "../base/geom_utils.h"



%ignore read_int;
%ignore read_double;

%include "../base/utils.h"


%include "../base/math_utils.h"
%include "../base/info.h"
%include "../base/bbox.h"
//%include "../base/timing.h"
%include "../base/polygons.h"
%include "../base/std_polys.h"

%include "../base/disp_poly.h"
%include "../base/vrml_writer.h"
%include "../base/pov_writer.h"
%include "../base/gl_writer.h"

