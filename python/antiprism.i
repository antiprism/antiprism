
%module antiprism

%{
#define SWIG_FILE_WITH_INIT
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

namespace std {
%template (VectorInt) vector<int>;
%template (VectorVectorInt) vector<vector<int> >;
%template (MapIntInt) map<int, int>;
%template (VectorVec3d) vector<anti::vec3d>;
}

%include "../base/antiprism.h"
%include "../base/const.h"
%include "../base/rand_gen.h"

%inline %{
   namespace anti {
      anti::col_geom_v read_geom(string file_name, char *err_msg)
         { anti::col_geom_v g; g.read(file_name, err_msg); return g; }
      anti::col_geom_v read_geom(char *err_msg)
         { return anti::read_geom("", err_msg); }
      anti::vec3d read_vec3d(char *str, char *err_msg)
         { anti::vec3d v; v.read(str, err_msg); return v; }
      anti::vec4d read_vec4d(char *str, char *err_msg)
         { anti::vec4d v; v.read(str, err_msg); return v; }
      anti::col_val read_color(char * str, char *err_msg)
         { anti::col_val col; col.read(str, err_msg); return col; }
   }
%}



%rename(__iadd__) anti::vec3d::operator+=;
%rename(__isub__) anti::vec3d::operator-=;
%rename(__imul__) anti::vec3d::operator*=;
%rename(__idiv__) anti::vec3d::operator/=;

%ignore anti::vec3d::read(char *);
%ignore anti::vec3d::operator[];
%ignore operator +(anti::vec3d, anti::vec3d);
%ignore operator -(anti::vec3d, anti::vec3d);
%ignore operator-(anti::vec3d);
%ignore operator -(anti::vec3d, anti::vec3d);
%ignore operator *(anti::vec3d, double);
%ignore operator *(double, anti::vec3d);
%ignore operator /(anti::vec3d, double);
%ignore operator /(double, anti::vec3d);

%extend anti::vec3d {
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



%rename(__neg__) anti::vec4d::operator-();
%rename(__iadd__) anti::vec4d::operator+=;
%rename(__isub__) anti::vec4d::operator-=;
%rename(__imul__) anti::vec4d::operator*=;
%rename(__idiv__) anti::vec4d::operator/=;

%ignore anti::vec4d::read(char *);
%ignore anti::vec4d::operator[];
%ignore operator +(anti::vec4d, anti::vec4d);
%ignore operator -(anti::vec4d, anti::vec4d);
%ignore operator-(anti::vec4d);
%ignore operator -(anti::vec4d, anti::vec4d);
%ignore operator *(anti::vec4d, double);
%ignore operator *(double, anti::vec4d);
%ignore operator /(anti::vec4d, double);
%ignore operator /(double, anti::vec4d);

%extend anti::vec4d {
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



%rename(__imul__) anti::mat3d::operator*=;

%ignore anti::mat3d::operator[];
%ignore operator *(const anti::mat3d &, const anti::vec3d &);
%ignore operator *(const anti::vec3d &, const anti::mat3d &);
%ignore operator *(const anti::mat3d &m1, const anti::mat3d &m2);
%ignore operator <(const anti::mat3d &m1, const anti::mat3d &m2);

%extend anti::mat3d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      vec3d operator *(const vec3d &vec) { return *$self * vec; }
      //vec3d operator *(const mat3d &mat) { return mat * (*$self); }
      mat3d operator *(const mat3d &mat) { return *$self * mat; }
      bool operator <(const mat3d &mat) { return *$self < mat; }
}

%include "../base/mat3d.h"



%ignore line_plane_intersect(anti::vec3d Q, anti::vec3d n, anti::vec3d P0, anti::vec3d P1, int *);

%include "../base/vec_utils.h"



%ignore mat4d::operator[];
%ignore operator *(const anti::mat4d &, const anti::vec4d &);
%ignore operator *(const anti::vec4d &, const anti::mat4d &);
%ignore operator *(const anti::mat4d &, const anti::mat4d &);
%ignore operator *(double, const anti::mat4d &);
%ignore operator *(const anti::mat4d &, double);
%ignore operator +(const anti::mat4d &, const anti::mat4d &);
%rename(__imul__) anti::mat4d::operator*=;
%rename(__iadd__) anti::mat4d::operator+=;

%extend anti::mat4d {
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



%ignore anti::col_val::operator[];
%ignore anti::col_val::read(char *);
%ignore anti::col_val::from_offvals(vector<char *> &);
%ignore anti::col_val::from_decvals(vector<double> &);
%ignore anti::col_val::from_intvals(vector<int> &);
%ignore anti::col_val::read_decvals(vector<char *> &);
%ignore anti::col_val::read_intvals(vector<char *> &);
%ignore anti::col_val::read_decvals(char *);
%ignore anti::col_val::read_intvals(char *);
%ignore anti::col_val::read_hexvals(char *);
%ignore anti::col_val::read_colorname(char *);

%extend anti::col_val {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
}

%include "../base/col_val.h"



%include "../base/col_geom.h"



%rename(is_set) anti::geom_if::operator bool() const;
%rename(__eq__) anti::col_geom_v::operator=(const anti::geom_if &);

%ignore anti::geom_if::add_hull(string);
%ignore anti::geom_if::add_hull(string);
%ignore anti::geom_if::set_hull(string);
%ignore anti::geom_if::read(string);
%ignore anti::geom_if::write(string);
%ignore anti::geom_if::write_crds(string);
%ignore anti::geom_if::read(FILE *);
%ignore anti::geom_if::read(FILE *, char *);
%ignore anti::geom_if::write(FILE *);
%ignore anti::geom_if::write(FILE *, char *);
%ignore anti::geom_if::write(FILE *, char *, int);

%extend anti::geom_if {
   //write() { return $self->write(""); }
   std::vector<anti::vec3d> allverts() { return $self->verts(); }
}

%include "../base/geom.h"

%extend anti::col_geom_v {
   void del_verts(std::vector<int> &v_idxs) { $self->delete_verts(v_idxs); }
}


%rename(__eq__) anti::color_map::operator=;

%ignore anti::color_map::read(string);

%include "../base/col_map.h"



%include "../base/coloring.h"



%rename(is_set) anti::t_set::operator bool() const;

%ignore operator *(const anti::t_set &, const anti::t_set &);
%ignore operator *=(anti::t_set &, const anti::t_set &);
%ignore operator +(const anti::t_set &, const anti::mat3d &);
%ignore operator +=(anti::t_set &, const anti::mat3d &);
%ignore anti::sch::read(string);

%extend anti::t_set {
   public:
      t_set operator *(const t_set &s) { return anti::t_set().product(*$self, s); }
      t_set &operator *=(const t_set &s) { return $self->product_with(s); }
      t_set operator +(const mat3d &m) { return anti::t_set(*$self).add(m); }
      t_set &operator +=(const mat3d &m) { return $self->add(m); }
}

anti::t_set &operator *=(t_set &s1, const t_set &s2);

%include "../base/symmetry.h"


%rename(__eq__) anti::scene_geom::operator=(const anti::scene_geom &);

%rename(__eq__) anti::scene::operator=(const anti::scene &);

%include "../base/normals.h"
%include "../base/bbox.h"
%include "../base/scene.h"


%include "../base/transforms.h"

%ignore anti::sort_merge_elems(anti::geom_if &geom, const string &merge_elems, double eps);


%ignore anti::minimum_distance(const anti::geom_if &geom, double);

%include "../base/geom_utils.h"



%ignore anti::read_int;
%ignore anti::read_double;

%include "../base/utils_ultragetopt.h"
%include "../base/utils.h"


%include "../base/math_utils.h"
%include "../base/info.h"
%include "../base/polygons.h"

%include "../base/disp_poly.h"
%include "../base/vrml_writer.h"
%include "../base/pov_writer.h"

