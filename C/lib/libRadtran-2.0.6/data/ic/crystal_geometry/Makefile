#
# Makefile for libRadtran's data/ic/crystal_geometry directory
#
# Bernhard Mayer
# bernhard.mayer@dlr.de

TOPDIR = ../../../

srcdir = .


THISDIR = crystal_geometry

DATAFILES = columns.crystal_points hexagon.crystal_points platesL40a40.crystal_points columnsL50a20.crystal_points plates.crystal_points

DISTFILES = Makefile.in $(DATAFILES)

MYSTICDISTFILES = $(DISTFILES)

include $(TOPDIR)/Makeconf

all: 

dist: all
	ln $(DISTFILES) ../../../`cat ../../../.fname`/data/ic/$(THISDIR)
.PHONY: dist

mystic-dist: all
	ln $(MYSTICDISTFILES) ../../../`cat ../../../.fname`/data/ic/$(THISDIR)
.PHONY: mystic-dist

install:
	$(MKPATH) $(datafiledir)/ic/$(THISDIR)
	$(INSTALL_DATA) $(DATAFILES) $(datafiledir)/ic/$(THISDIR)
.PHONY: clean

clean:

mostlyclean: clean

distclean: clean
	rm -f Makefile config.status

realclean: distclean

Makefile: Makefile.in ../../../config.status
	cd ../../../; sh config.status

# Prevent GNU make v3 from overflowing arg limit on SysV.
.NOEXPORT:

# eof
