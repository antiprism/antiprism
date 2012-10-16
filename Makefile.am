
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = base/tesselator base/expreval base src src_extra
if BUILD_ANTIVIEW
SUBDIRS += aview
endif

if BUILD_PY_MOD
# SUBDIRS += python
endif

EXTRA_DIST = bootstrap pad_file.xml README.git doc doc_src share

#             python/antiprism.i python/antiprism_doc.i \
#				 python/antiprism_wrap.cc python/antiprism.py \
#             python/distutils/antiprism_doc.i python/distutils/antiprism.py \
#             python/distutils/AUTHORS python/distutils/MANIFEST.in \
#             python/distutils/antiprism.i python/distutils/MANIFEST \
#             python/distutils/antiprism_wrap.cc python/distutils/COPYING
#
distcleancheck_listfiles = \
        find -type f -exec sh -c 'test -f $(srcdir)/{} || echo {}' ';'


supdir = $(datadir)/$(PACKAGE)-$(VERSION)
docdir = $(datadir)/doc/$(PACKAGE)-$(VERSION)
doc_DATA = README README.git NEWS AUTHORS COPYING

install-data-local:
	if test -d $(srcdir)/doc ; then \
	cd $(srcdir)/doc; \
	$(mkinstalldirs) $(DESTDIR)$(docdir); \
	for f in * css/* images/* lib_doc/html/* ; do \
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
	doxygen
	##cd python && doxygen && ./doxy2swig.py xml/index.xml antiprism_doc.i && rm -rf xml


