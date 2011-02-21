
LDADD = $(top_builddir)/base/libantiprism.la

#if USE_RPATH_FLAGS
#AM_LDFLAGS = -Wl,-rpath,$(libdir)
#endif


bin_PROGRAMS = antiview
antiview_SOURCES = antiview.cc vw_glut.cc \
                   vw_glut.h
antiview_LDADD = $(GLUT_LIB) $(LDADD)



