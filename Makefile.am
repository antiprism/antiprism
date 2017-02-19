
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = base/tesselator base/qhull base/muparser base src src_extra
if BUILD_ANTIVIEW
SUBDIRS += aview
endif

if BUILD_PY_MOD
SUBDIRS += python
endif

EXTRA_DIST = bootstrap pad_file.xml README.git Doxyfile vers_script.txt \
				 doc doc_src share debian \
             python/antiprism.i python/antiprism_doc.i \
				 python/antiprism_wrap.cc python/antiprism.py \
             python/distutils/antiprism_doc.i python/distutils/antiprism.py \
             python/distutils/AUTHORS python/distutils/MANIFEST.in \
             python/distutils/antiprism.i python/distutils/MANIFEST \
             python/distutils/antiprism_wrap.cc python/distutils/COPYING

distcleancheck_listfiles = \
        find -type f -exec sh -c 'test -f $(srcdir)/{} || echo {}' ';'


supdir = @datadir@/@PACKAGE@
docdir = @docdir@

doc_data_common = README README.git NEWS AUTHORS
if DEBIAN_INSTALL
doc_DATA = $(doc_data_common)
else
doc_DATA = $(doc_data_common) COPYING
endif

install-data-local:
	if test -d $(srcdir)/doc ; then \
	cd $(srcdir)/doc; \
	$(mkinstalldirs) $(DESTDIR)$(docdir); \
	for f in * css/* images/* ; do \
	## Find the install directory name
	d="`echo $$f | sed -e 's,[^/]*$$,,'`"; \
	## Make the install directory.
	test -d "$$d" && $(mkinstalldirs) $(DESTDIR)$(docdir)/$$d; \
	## Install the file
	test -f "$$f" && $(INSTALL_DATA) $$f $(DESTDIR)$(docdir)/$$f; \
	done \
	fi 
	cd $(srcdir)/share; \
	$(mkinstalldirs) $(DESTDIR)$(supdir); \
	for f in * col_lights/* col_maps/* extras/* extras/shell_scripts/* \
		extras/python/* models/* pov_inc/*; do \
	## Find the install directory name
	d="`echo $$f | sed -e 's,[^/]*$$,,'`"; \
	## Make the install directory.
	test -d "$$d" && $(mkinstalldirs) $(DESTDIR)$(supdir)/$$d; \
	## Install the file
	test -f "$$f" && $(INSTALL_DATA) $$f $(DESTDIR)$(supdir)/$$f; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)
	rm -rf $(DESTDIR)$(supdir)

prepare_release:
	cd doc_src && ./gtml doc.gtp && rm tmp.txt
	cd python && doxygen && ./doxy2swig.py xml/index.xml antiprism_doc.i && rm -rf xml

format_all:
	for f in base/*.cc base/*.h src/*.cc src/*.h aview/*.cc aview/*.h src_extra/*.cc ; do \
	clang-format -style=file -i $$f; \
	done

