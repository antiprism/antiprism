
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
%template (VectorVec3d) vector<anti::Vec3d>;
}

%include "../base/antiprism.h"
%include "../base/const.h"
%include "../base/random.h"

%inline %{
   namespace anti {
      anti::col_geom_v read_geom(string file_name, char *err_msg)
         { anti::col_geom_v g; g.read(file_name, err_msg); return g; }
      anti::col_geom_v read_geom(char *err_msg)
         { return anti::read_geom("", err_msg); }
      anti::Vec3d read_Vec3d(char *str, char *err_msg)
         { anti::Vec3d v; v.read(str, err_msg); return v; }
      anti::Vec4d read_Vec4d(char *str, char *err_msg)
         { anti::Vec4d v; v.read(str, err_msg); return v; }
      anti::Color read_color(char * str, char *err_msg)
         { anti::Color col; col.read(str, err_msg); return col; }
   }
%}



%rename(__iadd__) anti::Vec3d::operator+=;
%rename(__isub__) anti::Vec3d::operator-=;
%rename(__imul__) anti::Vec3d::operator*=;
%rename(__idiv__) anti::Vec3d::operator/=;

%ignore anti::Vec3d::read(char *);
%ignore anti::Vec3d::operator[];
%ignore operator +(anti::Vec3d, anti::Vec3d);
%ignore operator -(anti::Vec3d, anti::Vec3d);
%ignore operator-(anti::Vec3d);
%ignore operator -(anti::Vec3d, anti::Vec3d);
%ignore operator *(anti::Vec3d, double);
%ignore operator *(double, anti::Vec3d);
%ignore operator /(anti::Vec3d, double);
%ignore operator /(double, anti::Vec3d);

%extend anti::Vec3d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      Vec3d operator +(Vec3d v) { return *$self + v; }
      Vec3d operator -(Vec3d v) { return *$self - v; }
      Vec3d operator -() { return -(*$self); }
      Vec3d operator *(double n) { return *$self * n; }
      Vec3d operator /(double n) { return *$self / n; }
}

%include "../base/vec3d.h"



%rename(__neg__) anti::Vec4d::operator-();
%rename(__iadd__) anti::Vec4d::operator+=;
%rename(__isub__) anti::Vec4d::operator-=;
%rename(__imul__) anti::Vec4d::operator*=;
%rename(__idiv__) anti::Vec4d::operator/=;

%ignore anti::Vec4d::read(char *);
%ignore anti::Vec4d::operator[];
%ignore operator +(anti::Vec4d, anti::Vec4d);
%ignore operator -(anti::Vec4d, anti::Vec4d);
%ignore operator-(anti::Vec4d);
%ignore operator -(anti::Vec4d, anti::Vec4d);
%ignore operator *(anti::Vec4d, double);
%ignore operator *(double, anti::Vec4d);
%ignore operator /(anti::Vec4d, double);
%ignore operator /(double, anti::Vec4d);

%extend anti::Vec4d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      Vec4d operator +(Vec4d v) { return *$self + v; }
      Vec4d operator -(Vec4d v) { return *$self - v; }
      Vec4d operator -() { return -(*$self); }
      Vec4d operator *(double n) { return *$self * n; }
      Vec4d operator /(double n) { return *$self / n; }
}

%include "../base/vec4d.h"



%rename(__imul__) anti::Trans3d::operator*=;

%ignore anti::Trans3d::operator[];
%ignore operator *(const anti::Trans3d &, const anti::Vec3d &);
%ignore operator *(const anti::Vec3d &, const anti::Trans3d &);
%ignore operator *(const anti::Trans3d &m1, const anti::mat3d &m2);
%ignore operator <(const anti::Trans3d &m1, const anti::mat3d &m2);

%extend anti::Trans3d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      Vec3d operator *(const Vec3d &vec) { return *$self * vec; }
      //Vec3d operator *(const Trans3d &mat) { return mat * (*$self); }
      Trans3d operator *(const mat3d &mat) { return *$self * mat; }
      bool operator <(const Trans3d &mat) { return *$self < mat; }
}

%include "../base/trans3d.h"



%ignore line_plane_intersect(anti::Vec3d Q, anti::Vec3d n, anti::Vec3d P0, anti::Vec3d P1, int *);

%include "../base/vec_utils.h"



%ignore Trans4d::operator[];
%ignore operator *(const anti::Trans4d &, const anti::Vec4d &);
%ignore operator *(const anti::Vec4d &, const anti::Trans4d &);
%ignore operator *(const anti::Trans4d &, const anti::Trans4d &);
%ignore operator *(double, const anti::Trans4d &);
%ignore operator *(const anti::Trans4d &, double);
%ignore operator +(const anti::Trans4d &, const anti::Trans4d &);
%rename(__imul__) anti::Trans4d::operator*=;
%rename(__iadd__) anti::Trans4d::operator+=;

%extend anti::Trans4d {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
      void __setitem__(int idx, double val) { $self->operator[](idx) = val; }
      Vec4d operator *(const Vec4d &vec) { return *$self * vec; }
      //Vec4d operator *(const Trans4d &mat) { return *$self * mat; }
      Trans4d operator *(const Trans4d &mat) { return *$self * mat; }
      Trans4d operator *(double n) { return *$self * n; }
      Trans4d operator +(const Trans4d &mat) { return *$self + mat; }
}

%include "../base/trans4d.h"



%ignore anti::Color::operator[];
%ignore anti::Color::read(char *);
%ignore anti::Color::from_offvals(vector<char *> &);
%ignore anti::Color::from_decvals(vector<double> &);
%ignore anti::Color::from_intvals(vector<int> &);
%ignore anti::Color::read_decvals(vector<char *> &);
%ignore anti::Color::read_intvals(vector<char *> &);
%ignore anti::Color::read_decvals(char *);
%ignore anti::Color::read_intvals(char *);
%ignore anti::Color::read_hexvals(char *);
%ignore anti::Color::read_colorname(char *);

%extend anti::Color {
   public:
      double __getitem__(int idx) { return $self->operator[](idx); }
}

%include "../base/color.h"



%rename(is_set) anti::Geometry::operator bool() const;
%rename(__eq__) anti::Geometry::operator=(const anti::Geometry &);

%ignore anti::Geometry::add_hull(string);
%ignore anti::Geometry::add_hull(string);
%ignore anti::Geometry::set_hull(string);
%ignore anti::Geometry::read(string);
%ignore anti::Geometry::write(string);
%ignore anti::Geometry::write_crds(string);
%ignore anti::Geometry::read(FILE *);
%ignore anti::Geometry::read(FILE *, char *);
%ignore anti::Geometry::write(FILE *);
%ignore anti::Geometry::write(FILE *, char *);
%ignore anti::Geometry::write(FILE *, char *, int);

%extend anti::Geometry {
   //write() { return $self->write(""); }
   std::vector<anti::Vec3d> allverts() { return $self->verts(); }
}

%include "../base/elemprops.h"
%include "../base/geometry.h"

%extend anti::Geometry {
   void del_verts(std::vector<int> &v_idxs) { $self->delete_verts(v_idxs); }
}


%rename(__eq__) anti::color_map::operator=;

%ignore anti::ColorMap::read(string);

%include "../base/colormap.h"



%include "../base/coloring.h"



%rename(is_set) anti::t_set::operator bool() const;

%ignore operator *(const anti::t_set &, const anti::t_set &);
%ignore operator *=(anti::t_set &, const anti::t_set &);
%ignore operator +(const anti::t_set &, const anti::Trans3d &);
%ignore operator +=(anti::t_set &, const anti::Trans3d &);
%ignore anti::sch::read(string);

%extend anti::t_set {
   public:
      t_set operator *(const t_set &s) { return anti::t_set().product(*$self, s); }
      t_set &operator *=(const t_set &s) { return $self->product_with(s); }
      t_set operator +(const Trans3d &m) { return anti::t_set(*$self).add(m); }
      t_set &operator +=(const Trans3d &m) { return $self->add(m); }
}

anti::t_set &operator *=(t_set &s1, const t_set &s2);

%include "../base/symmetry.h"


%rename(__eq__) anti::scene_geom::operator=(const anti::scene_geom &);

%rename(__eq__) anti::scene::operator=(const anti::scene &);

%include "../base/normal.h"
%include "../base/boundbox.h"
%include "../base/scene.h"



%ignore anti::sort_merge_elems(anti::Geometry &geom, const string &merge_elems, double eps);


%ignore anti::minimum_distance(const anti::Geometry &geom, double);

%include "../base/geometryutils.h"



%ignore anti::read_int;
%ignore anti::read_double;

%include "../base/getopt.h"
%include "../base/utils.h"


%include "../base/mathutils.h"
%include "../base/geometryinfo.h"
%include "../base/polygon.h"

%include "../base/programopts.h"
%include "../base/displaypoly.h"
%include "../base/vrmlwriter.h"
%include "../base/povwriter.h"

