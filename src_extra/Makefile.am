
#LDADD = ../base/libantiprism.la
LDFLAGS += -L../base -lantiprism

if USE_RPATH_FLAGS
AM_LDFLAGS = -Wl,-rpath,$(libdir)
endif


bin_PROGRAMS = twist spidron off2txt dome_layer

twist_SOURCES = twist.cc
spidron_SOURCES = spidron.cc
off2txt_SOURCES = off2txt.cc
dome_layer_SOURCES = dome_layer.cc

