
LDADD = $(top_builddir)/base/libantiprism.la
#LDFLAGS += -L../base -lantiprism

#if USE_RPATH_FLAGS
#AM_LDFLAGS = -Wl,-rpath,$(libdir)
#endif

.cc.1: $(top_srcdir)/configure.ac
	$(MAKE) $(AM_MAKEFLAGS) $*$(EXEEXT)
	grep "doc/$*" ../doc_src/doc.gtp | sed -e 's/\(.*\) 2 \(.*\)/[NAME]\n\2/' > /tmp/h2m_name_$*
	help2man -i /tmp/h2m_name_$* --output=$@ ./$*$(EXEEXT)
	rm /tmp/h2m_name_$*


noinst_LTLIBRARIES = liblattice_grid.la
					 
liblattice_grid_la_SOURCES = lattice_grid.cc lattice_grid.h

bin_PROGRAMS = off2pov off2m off2vrml off2crds m2off \
      off_color off_util off_trans off_align \
      poly_kscope polygon zono conv_hull pol_recip \
		geodesic minmax sph_saff sph_rings \
		poly_weave off_report off_query \
		kcycle string_art jitterbug unitile2d repel \
		lat_grid lat_util canonical conway n_icons iso_delta \
		bravais waterman col_util planar

dist_man1_MANS = off2pov.1 off2m.1 off2vrml.1 off2crds.1 m2off.1 \
      off_color.1 off_util.1 off_trans.1 off_align.1 \
      poly_kscope.1 polygon.1 zono.1 conv_hull.1 pol_recip.1 \
		geodesic.1 minmax.1 sph_saff.1 sph_rings.1 \
		poly_weave.1 off_report.1 off_query.1 \
		kcycle.1 string_art.1 jitterbug.1 unitile2d.1 repel.1 \
		lat_grid.1 lat_util.1 canonical.1 conway.1 n_icons.1 iso_delta.1 \
		bravais.1 waterman.1 col_util.1 planar.1

off2pov_SOURCES = off2pov.cc
off2crds_SOURCES = off2crds.cc
off2vrml_SOURCES = off2vrml.cc
m2off_SOURCES = m2off.cc
off2m_SOURCES = off2m.cc
off_color_SOURCES = off_color.cc
off_util_SOURCES = off_util.cc help.h
off_trans_SOURCES = off_trans.cc
off_align_SOURCES = off_align.cc
poly_kscope_SOURCES = poly_kscope.cc
polygon_SOURCES = polygon.cc
zono_SOURCES = zono.cc
conv_hull_SOURCES = conv_hull.cc
pol_recip_SOURCES = pol_recip.cc
geodesic_SOURCES = geodesic.cc
minmax_SOURCES = minmax.cc
sph_saff_SOURCES = sph_saff.cc
sph_rings_SOURCES = sph_rings.cc
poly_weave_SOURCES = poly_weave.cc
off_report_SOURCES = off_report.cc rep_print.cc rep_print.h
off_query_SOURCES = off_query.cc rep_print.cc rep_print.h
kcycle_SOURCES = kcycle.cc
string_art_SOURCES = string_art.cc
jitterbug_SOURCES = jitterbug.cc
unitile2d_SOURCES = unitile2d.cc
repel_SOURCES = repel.cc
lat_grid_SOURCES = lat_grid.cc
lat_grid_LDADD = liblattice_grid.la $(LDADD)
lat_util_SOURCES = lat_util.cc
lat_util_LDADD = liblattice_grid.la $(LDADD)
canonical_SOURCES = canonical.cc
conway_SOURCES = conway.cc
conway_LDADD = liblattice_grid.la $(LDADD)
n_icons_SOURCES = n_icons.cc n_icons.h
iso_delta_SOURCES = iso_delta.cc
bravais_SOURCES = bravais.cc
bravais_LDADD = liblattice_grid.la $(LDADD)
waterman_SOURCES = waterman.cc
waterman_LDADD = liblattice_grid.la $(LDADD)
col_util_SOURCES = col_util.cc
planar_SOURCES = planar.cc

