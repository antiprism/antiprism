
BUILT_SOURCES = antiprism_wrap.cc

if USE_SWIG
SWIG_SRCS = antiprism.i antiprism_doc.i
SWIG_CMD = $(SWIG) -c++ $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/base \
                      -o antiprism_wrap.cc $(top_srcdir)/python/antiprism.i
else
  SWIG_CMD = echo "swig is not installed so antiprism_wrap.cc is not regenerated" ; cp $(top_srcdir)/python/antiprism_wrap.cc tmp.cc ; cp tmp.cc antiprism_wrap.cc;
endif


#pkgpython_PYTHON = antiprism.py
#pkgpyexec_LTLIBRARIES = _antiprism.la
#_antiprism_la_SOURCES = antiprism_wrap.cc $(SWIG_SRCS)
#_antiprism_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/base
#_antiprism_la_LDFLAGS = -module
#_antiprism_la_LIBADD = ../base/libantiprism.la


../base/libantiprism.la:

#	$(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/base -o antiprism_wrap.cc \
#	                                      $(top_srcdir)/python/antiprism.i

antiprism_wrap.cc : $(SWIG_SRCS) ../base/libantiprism.la
	$(SWIG_CMD)
	rm -rf distutils/lib_src \
	       distutils/build \
	       distutils/dist
	mkdir distutils/lib_src \
	      distutils/lib_src/muparser \
			distutils/lib_src/qhull \
			distutils/lib_src/tesselator
	cp $(top_srcdir)/base/*.h $(top_srcdir)/base/*.cc \
	                                   distutils/lib_src
	cp $(top_srcdir)/base/muparser/*.h $(top_srcdir)/base/muparser/*.cc \
	                                   distutils/lib_src/muparser
	cp $(top_srcdir)/base/qhull/*.h $(top_srcdir)/base/qhull/*.cc \
	                                   distutils/lib_src/qhull
	cp $(top_srcdir)/base/tesselator/*.h $(top_srcdir)/base/tesselator/*.cc \
	                                   distutils/lib_src/tesselator
	sed -e 's/\.\.\/base/lib_src/' $(top_srcdir)/python/antiprism.i > \
	                                   distutils/antiprism.i
	sed -e 's/\.\.\/base/lib_src/' $(top_srcdir)/python/antiprism_wrap.cc > \
	                                   distutils/antiprism_wrap.cc
	cp $(top_srcdir)/python/antiprism.py     distutils
	cp $(top_srcdir)/python/antiprism_doc.i  distutils
	cp $(top_srcdir)/AUTHORS                 distutils
	cp $(top_srcdir)/COPYING                 distutils

clean-local:
	rm -rf distutils/build \
	       distutils/dist \
	       distutils/lib_src \
	       distutils/antiprism.i \
			 distutils/antiprism_doc.i \
	       distutils/antiprism_wrap.cc \
			 distutils/antiprism.py \
	       distutils/AUTHORS \
			 distutils/COPYING \
			 distutils/MANIFEST \
	       distutils/antiprism_wrap.cc \
			 distutils/antiprism.py
	if test -f tmp.cc ; then  rm -f antiprism_wrap.cc ; fi
	rm -f tmp.cc
	




	
