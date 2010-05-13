
lib_LTLIBRARIES = libantiprism.la

libantiprism_la_SOURCES = off_read.cc off_write.cc crds_read.cc disp_poly.cc\
   geom.cc col_geom.cc col_map.cc col_val.cc dual.cc vec3d.cc mat3d.cc \
   vec4d.cc mat4d.cc vec_utils.cc vec_utils_norm.cc vec_utils_cent.cc \
	utils.cc math_utils.cc c_hull.cc triangulate.cc symmetry.cc \
	sort_merge.cc bbox.cc col_read.cc info.cc \
   coloring.cc prop_col.cc named_cols.cc geodesic.cc zonohedron.cc \
	johnson.cc uniform.cc kaleido.cc std_polys.cc skilling.cc \
	timing.cc polygons.cc pov_writer.cc scene.cc \
	gl_writer.cc stipple.cc gl_stroke_font.cc \
	canonic.cc trans.cc faces.cc vrml_writer.cc \
   \
   qh_user.cc qh_global.cc qh_stat.cc qh_io.cc \
   qh_geom2.cc qh_poly2.cc qh_merge.cc qh_qhull.cc \
	qh_geom.cc qh_poly.cc qh_qset.cc qh_mem.cc \
   \
	antiprism.h bbox.h col_map.h const.h disp_poly.h\
	geodesic.h geom.h geom_utils.h gl_stroke_font.h gl_writer.h \
	col_val.h col_geom.h info.h mat3d.h mat4d.h math_utils.h named_cols.h \
   off_file.h pov_writer.h coloring.h prop_col.h timing.h \
	polygons.h std_polys.h skilling.h kaleido.h \
   \
	qh_geom.h qh_io.h qh_mem.h qh_merge.h qh_poly.h qh_qhull_a.h \
	qh_qhull.h qh_qset.h qh_stat.h qh_user.h \
	\
	rand_gen.h scene.h symmetry.h transforms.h \
	utils.h vec3d.h vec4d.h vec_utils.h vrml_writer.h

supdir = $(datadir)/$(PACKAGE)-$(VERSION)
libantiprism_la_CPPFLAGS = -DSUPDIR="\"$(supdir)\""
libantiprism_la_LDFLAGS = -no-undefined -version-info 0:0:0

pkginclude_HEADERS = antiprism.h \
						const.h vec3d.h mat3d.h vec_utils.h vec4d.h \
						mat4d.h col_val.h col_map.h col_geom.h geom.h coloring.h \
						symmetry.h transforms.h geom_utils.h math_utils.h utils.h \
						info.h bbox.h timing.h polygons.h std_polys.h skilling.h kaleido.h \
						scene.h pov_writer.h vrml_writer.h gl_writer.h


