
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = base/tesselator base src src_extra
if BUILD_ANTIVIEW
SUBDIRS += aview
endif

if BUILD_PY_MOD
# SUBDIRS += python
endif

EXTRA_DIST = bootstrap pad_file.xml doc share base/qh_COPYING.txt 

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
doc_DATA = README NEWS AUTHORS COPYING

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
	mkdir -p doc/man
	cd doc_src && ./gtml doc.gtp && rm tmp.txt
	doxygen
	##cd python && doxygen && ./doxy2swig.py xml/index.xml antiprism_doc.i && rm -rf xml


man1_MANS = doc/man/antiview.1 doc/man/bravais.1 doc/man/canonical.1 \
doc/man/col_util.1 doc/man/conv_hull.1 doc/man/conway.1 \
doc/man/dome_layer.1 doc/man/geodesic.1 doc/man/iso_delta.1 \
doc/man/jitterbug.1 doc/man/kcycle.1 doc/man/lat_grid.1 \
doc/man/lat_util.1 doc/man/m2off.1 doc/man/minmax.1 \
doc/man/n_icons.1 doc/man/off2crds.1 doc/man/off2m.1 \
doc/man/off2pov.1 doc/man/off2txt.1 doc/man/off2vrml.1 \
doc/man/off_align.1 doc/man/off_color.1 doc/man/off_normals.1 \
doc/man/off_query.1 doc/man/off_report.1 doc/man/off_trans.1 \
doc/man/off_util.1 doc/man/planar.1 doc/man/pol_recip.1 \
doc/man/polygon.1 doc/man/poly_kscope.1 doc/man/poly_weave.1 \
doc/man/repel.1 doc/man/sph_rings.1 doc/man/sph_saff.1 \
doc/man/spidron.1 doc/man/string_art.1 doc/man/twist.1 \
doc/man/unitile2d.1 doc/man/waterman.1 doc/man/zono.1


