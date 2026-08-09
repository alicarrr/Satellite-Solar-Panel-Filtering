#
# Makefile for libRadtran's example/triangular_srfc directory
#

TOPDIR = ../../..

srcdir = .


THISDIR = $(notdir $(CURDIR))

SUBDIRS =

EXAMPLEFILES = \
	atm.dat           \
  triangle_input.nc \
	no_atmosphere.inp \
	atmosphere.inp \
	99_flat_surface_thermal_atmosphere.testuite_output.nc \
	99_flat_surface_thermal_no_atmosphere.testuite_output.nc \

DISTFILES = Makefile.in $(EXAMPLEFILES)

include $(TOPDIR)/Makeconf

all:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: all

dist: all
	for dir in $(SUBDIRS); do mkdir $(TOPDIR)/`cat $(TOPDIR)/.fname`/examples/triangular_srfc/$(THISDIR)/$$dir; $(MAKE) -C $$dir $@ || exit $$?; done; \
	ln $(DISTFILES) $(TOPDIR)/`cat $(TOPDIR)/.fname`/examples/triangular_srfc/$(THISDIR)
.PHONY: dist

mystic-dist: dist
.PHONY: mystic-dist

install: all
	$(MKPATH) $(examplesdir)/$(THISDIR)
	$(INSTALL_DATA) $(EXAMPLEFILES) $(examplesdir)/$(THISDIR)
.PHONY: install

clean:
	rm -f mc* randomseed *.xmf hdf5_compatible_mc.flx.triangle.nc
.PHONY: clean

mostlyclean: clean
.PHONY: mostlyclean

distclean: mostlyclean
	rm -f Makefile
.PHONY: distclean

realclean: distclean
.PHONY: realclean

Makefile: Makefile.in $(TOPDIR)/config.status
	cd $(TOPDIR); sh config.status
# end default targets

uvspec_bin:
	make -C $(TOPDIR) uvspec
.PHONY: uvspec_bin

triangle_input.nc: ../02_gen_flat_surface.py
	./$< $@ -half_area -srfc_height 1e-6 -T 300

triangle_input.nc.xmf: triangle_input.nc
	python ${TOPDIR}/external/xmftools/xmftool.py $<

mc.withatm.flx.triangle.nc: atmosphere.inp uvspec_bin
	make -C $(TOPDIR) uvspec
	$(TOPDIR)/bin/uvspec -f $<

mc.noatm.flx.triangle.nc: no_atmosphere.inp uvspec_bin
	make -C $(TOPDIR) uvspec
	$(TOPDIR)/bin/uvspec -f $<
