if DEBIAN_INSTALL
pkglib_LTLIBRARIES = libantiprism.la
else
lib_LTLIBRARIES = libantiprism.la
endif

libantiprism_la_SOURCES = \
	off_read.cc off_write.cc crds_read.cc displaypoly.cc\
	geometry.cc geometryutils.cc colormap.cc color.cc dual.cc \
	programopts.cc status.cc vec3d.cc trans3d.cc \
	vec4d.cc trans4d.cc vec_utils.cc vec_utils_norm.cc vec_utils_cent.cc \
	utils.cc utils_parser.cc getopt.cc mathutils.cc \
	normal.cc c_hull.cc triangulate.cc \
	symmetry.cc sort_merge.cc boundbox.cc geometryinfo.cc \
	coloring.cc prop_col.cc named_cols.cc geodesic.cc zonohedron.cc \
	johnson.cc uniform.cc std_polys.cc skilling.cc \
	timer.cc polygon.cc povwriter.cc scene.cc \
	canonic.cc trans.cc faces.cc vrmlwriter.cc wythoff.cc \
	\
	antiprism.h boundbox.h elemprops.h colormap.h coloring.h color.h \
	const.h displaypoly.h geometry.h geometryutils.h geometryinfo.h \
	trans3d.h trans4d.h mathutils.h normal.h polygon.h povwriter.h \
	programopts.h random.h scene.h status.h symmetry.h timer.h utils.h \
	getopt.h vec3d.h vec4d.h vec_utils.h vrmlwriter.h \
	\
	private_geodesic.h private_misc.h private_named_cols.h \
	private_off_file.h private_prop_col.h private_std_polys.h

supdir = $(datadir)/$(PACKAGE)
libantiprism_la_CPPFLAGS = -DSUPDIR="\"$(supdir)\"" 
libantiprism_la_LDFLAGS = -no-undefined -version-info 0:0:0
if HAVE_LD_VERSION_SCRIPT
libantiprism_la_LDFLAGS += \
   -Wl,--version-script,'$(top_srcdir)/vers_script.txt'
endif
libantiprism_la_LIBADD = tesselator/libtess.la qhull/libqhull.la \
	muparser/libmuparser.la

if DEBIAN_INSTALL
pkginclude_HEADERS =
else
pkginclude_HEADERS = \
	antiprism.h \
	boundbox.h \
	colormap.h \
	color.h \
	coloring.h \
	const.h \
	displaypoly.h \
	geometry.h \
	geometryutils.h \
	geometryinfo.h \
	mathutils.h \
	normal.h \
	polygon.h \
	povwriter.h \
	programopts.h \
	elemprops.h \
	random.h \
	scene.h \
	status.h \
	symmetry.h \
	timer.h \
	trans3d.h \
	trans4d.h \
	utils.h \
	getopt.h \
	vec3d.h \
	vec4d.h \
	vec_utils.h \
	vrmlwriter.h
	
endif
