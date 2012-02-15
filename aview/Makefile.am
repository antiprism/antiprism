
LIBS = $(GLPROG_LIBS)
CXXFLAGS = $(GLPROG_CXXFLAGS)

.cc.1: $(top_srcdir)/configure.ac
	$(MAKE) $(AM_MAKEFLAGS) $*$(EXEEXT)
	help2man --output=$@ ./$*$(EXEEXT)


LDADD = libantiprism_dispgl.la $(top_builddir)/base/libantiprism.la

noinst_LTLIBRARIES = libantiprism_dispgl.la
libantiprism_dispgl_la_SOURCES=gl_writer.cc stipple.cc gl_stroke_font.cc \
	disp_poly_gl.cc \
	gl_stroke_font.h gl_writer.h disp_poly_gl.h

supdir = $(datadir)/$(PACKAGE)-$(VERSION)
libantiprism_dispgl_la_CPPFLAGS = -DSUPDIR="\"$(supdir)\""
libantiprism_dispgl_la_LDFLAGS = -no-undefined -version-info 0:0:0


bin_PROGRAMS = antiview

man1_MANS = antiview.1

antiview_SOURCES = antiview.cc vw_glut.cc \
                   vw_glut.h
antiview_LDADD = $(GLUT_LIBS) $(LDADD)
antiview_CXXFLAGS = $(GLUT_CFLAGS)



