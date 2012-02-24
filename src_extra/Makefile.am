
LDADD = $(top_builddir)/base/libantiprism.la
#LDFLAGS += -L../base -lantiprism

#if USE_RPATH_FLAGS
#AM_LDFLAGS = -Wl,-rpath,$(libdir)
#endif

.cc.1: $(top_srcdir)/configure.ac
	$(MAKE) $(AM_MAKEFLAGS) $*$(EXEEXT)
	grep "doc/$*" ../doc_src/doc.gtp | sed -e 's/\(.*\)UNSUPPORTED: \(.*\)/[NAME]\n\2/' > /tmp/h2m_name_$*
	help2man -i /tmp/h2m_name_$* --output=$@ ./$*$(EXEEXT)
	rm /tmp/h2m_name_$*


bin_PROGRAMS = twist spidron off2txt dome_layer off_normals

dist_man1_MANS = twist.1 spidron.1 off2txt.1 dome_layer.1 off_normals.1

twist_SOURCES = twist.cc
spidron_SOURCES = spidron.cc
off2txt_SOURCES = off2txt.cc
dome_layer_SOURCES = dome_layer.cc
off_normals_SOURCES = off_normals.cc

