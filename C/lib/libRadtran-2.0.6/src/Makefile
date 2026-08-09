#
# Makefile for libRadtran's src directory
#
# Arve Kylling
# arve.kylling@gmail.no

TOPDIR = ..

include $(TOPDIR)/Makeconf

%.c: %.l
	$(LEX) $(LFLAGS) -o $@ $<

srcdir = .


PYTHONDIR = $(TOPDIR)/src_py/

ADD_WORLDLOOP_RESULT = add_worldloop_result.c
ANGRES               = angres.c
SNOWALBEDO           = snowalbedo.c
VRSALBEDO            = vrs_ocean_albedo.c
CONVOLUTION          = conv.c
INTEGRATE            = integrate.c
SPLINE_INTERPOLATION = spline.c
FULIOU_FSOURCES      = ckdfuf.f ckdrayf.f ckdcldf.f ckdicef.f \
                       fl_dkfrwn_0599.f fl_misc_subs.f        \
                       fl_trmm_wnflt_0599.f fl_water_hu2.f
FULIOU_CSOURCES      = ckdfu.c
FULIOU_INCLUDES      = fl_radparams.cinc fl_radparams.inc \
                       fl_special.inc fl_trmm_window.inc
CLDPRP               = cldprp.c atmosphere.c cloud.c redistribute.c \
                       $(FULIOU_CSOURCES) ancillary.c  options_stash.c \
                       aerosol.c ck.c ipa.c
CLDGEN               = $(CLDGEN_TSOURCES)
PLKAVG               = plkavg.c
CLDPRP_DIST          = cldprp.c
TIME2SZA             = time2sza.c
SZA2TIME             = sza2time.c
ZENITH               = zenith.c
NOON                 = noon.c
MAKE_SLITFUNCTION    = make_slitfunction.c
MAKE_ANGRESFUNC      = make_angresfunc.c
DER_STAM             = read_Stamnes_tab.c
MIE                  = mie.c
PHASE                = phase.c
PMOM                 = pmom.c
MIE_HEADER           = $(patsubst %.c, %.h, $(MIE))
MIE_LEXSRC           = mie_lex.l
SOFI                 = sofi.c
ARLEM                = $(ARLEM_SOURCES)
UVSPEC               = uvspec.c ancillary.c options_stash.c atmosphere.c extraterrestrial.c \
                       aerosol.c cloud.c cloud3d.c elevation2d.c temperature2d.c \
                       molecular.c sza.c albedo.c \
                       redistribute.c fluorescence.c molecular3d.c $(TIPA_SOURCES) \
                       solve_rte.c ipa.c ck.c $(FULIOU_CSOURCES) $(CONDOR_SOURCES)
UVSPEC_FUNCTION      = uvspecfunction.c
UVSPEC_MPI           = uvspec_mpi.c
UVSPEC_FORTRAN       = $(FULIOU_FSOURCES) $(LOWTRAN_FSOURCES)
UVSPEC_HEADER        = $(patsubst %.c, %.h, $(UVSPEC))
UVSPEC_LEXSRC        = $(PYTHONDIR)uvspec_lex.l
VAC2AIR              = vac2air.c
SSRADAR              = ssradar.c

PERL_MODULES         = UVTools.pm

ADD_WORLDLOOP_RESULT_OBJECTS = $(patsubst %.c, %.o, $(ADD_WORLDLOOP_RESULT))
ANGRES_OBJECTS               = $(patsubst %.c, %.o, $(ANGRES))
SNOWALBEDO_OBJECTS           = $(patsubst %.c, %.o, $(SNOWALBEDO))
VRSALBEDO_OBJECTS            = $(patsubst %.c, %.o, $(VRSALBEDO))
CONVOLUTION_OBJECTS          = $(patsubst %.c, %.o, $(CONVOLUTION))
INTEGRATE_OBJECTS            = $(patsubst %.c, %.o, $(INTEGRATE))
SPLINE_INTERPOLATION_OBJECTS = $(patsubst %.c, %.o, $(SPLINE_INTERPOLATION))
CLDPRP_OBJECTS               = $(patsubst %.c, %.o, $(CLDPRP)) \
                               $(patsubst %.f, %.o, $(UVSPEC_FORTRAN))
CLDGEN_OBJECTS               = $(patsubst %.c, %.o, $(CLDGEN))
PLKAVG_OBJECTS               = $(patsubst %.c, %.o, $(PLKAVG))
TIME2SZA_OBJECTS             = $(patsubst %.c, %.o, $(TIME2SZA))
SZA2TIME_OBJECTS             = $(patsubst %.c, %.o, $(SZA2TIME))
ZENITH_OBJECTS               = $(patsubst %.c, %.o, $(ZENITH))
NOON_OBJECTS                 = $(patsubst %.c, %.o, $(NOON))
MAKE_SLITFUNCTION_OBJECTS    = $(patsubst %.c, %.o, $(MAKE_SLITFUNCTION))
MAKE_ANGRESFUNC_OBJECTS      = $(patsubst %.c, %.o, $(MAKE_ANGRESFUNC))
DER_STAM_OBJECTS             = $(patsubst %.c, %.o, $(DER_STAM))
MIE_OBJECTS                  = $(patsubst %.c, %.o, $(MIE)) $(patsubst %.l, %.o, $(MIE_LEXSRC))
PHASE_OBJECTS                = $(patsubst %.c, %.o, $(PHASE))
PMOM_OBJECTS                 = $(patsubst %.c, %.o, $(PMOM))
SOFI_OBJECTS                 = $(patsubst %.c, %.o, $(SOFI))
UVSPEC_OBJECTS               = $(patsubst %.c, %.o, $(UVSPEC))         \
                               $(patsubst %.f, %.o, $(UVSPEC_FORTRAN))
ARLEM_OBJECTS                = $(patsubst %.c, %.o, $(ARLEM))
UVSPEC_FUNCTION_OBJECTS      = $(patsubst %.c, %.o, $(UVSPEC_FUNCTION))
UVSPEC_MPI_OBJECTS           = $(patsubst %.c, %.o, $(UVSPEC_MPI))
VAC2AIR_OBJECTS              = $(patsubst %.c, %.o, $(VAC2AIR))
SSRADAR_OBJECTS              = $(patsubst %.c, %.o, $(SSRADAR))

RADTRAN_F_LIB =  ../lib/libRadtran_f.a
RADTRAN_C_LIB =  ../lib/libRadtran_c.a

# duplicate for the linker to resolve circular dependencies
RADTRAN_LIB = $(RADTRAN_C_LIB) $(RADTRAN_F_LIB) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)

DISTFILES = Makefile.in .depend \
	$(INTEGRATE) $(ANGRES) $(SNOWALBEDO) $(CONVOLUTION) \
	$(SPLINE_INTERPOLATION) $(CLDPRP_DIST) $(PLKAVG) \
        $(SZA2TIME) $(ZENITH) $(NOON) $(VRSALBEDO) \
	$(TIME2SZA) $(MAKE_SLITFUNCTION) $(MAKE_ANGRESFUNC) $(DER_STAM) \
	$(MIE) $(PHASE) $(PMOM) $(SOFI) $(UVSPEC) $(VAC2AIR) \
	$(MIE_LEXSRC) mie_lex.c $(MIE_HEADER)\
	$(UVSPEC_LEXSRC) uvspec_lex.c $(UVSPEC_HEADER) \
	$(FULIOU_FSOURCES) $(FULIOU_INCLUDES) \
	$(LOWTRAN_FSOURCES) $(LOWTRAN_PARAM)  \
	UVTools.pm Gen_o3_tab.pl Gen_wc_tab.pl Gen_snow_tab.pl MIE_SMOOTH \
        ndiff Calc_J.pl cloudprp2cdf.sh cloudprp2cdf.awk addlevel addlevel.c \
        airmass.pl \
        $(UVSPEC_FUNCTION) $(UVSPEC_MPI) $(SSRADAR)

MYSTICDISTFILES = $(DISTFILES)

SRC_C = $(ANGRES) $(SNOWALBEDO) $(CONVOLUTION) $(INTEGRATE) \
        $(SPLINE_INTERPOLATION) $(CLDPRP) $(CLDGEN) $(PLKAVG) $(CLDPRP_DIST) \
        $(TIME2SZA) $(SZA2TIME) $(ZENITH) $(NOON) $(MAKE_SLITFUNCTION) \
        $(MAKE_ANGRESFUNC) $(DER_STAM) $(MIE) mie_lex.c $(PHASE) $(PMOM) $(SOFI) \
	$(UVSPEC) $(ARLEM) $(UVSPEC_FUNCTION) $(UVSPEC_MPI) uvspec_lex.c $(VAC2AIR) \
	$(VRSALBEDO) $(SSRADAR)

BINARIES = \
					 $(ARLEM_TOOL) \
					 ../bin/angres \
					 ../bin/cldprp \
					 ../bin/conv \
					 ../bin/integrate \
					 ../bin/make_angresfunc \
					 ../bin/make_slitfunction \
					 ../bin/mie \
					 ../bin/noon \
					 ../bin/phase \
					 ../bin/plkavg \
					 ../bin/pmom \
					 ../bin/read_Stamnes_tab \
					 ../bin/snowalbedo \
					 ../bin/sofi \
					 ../bin/spline \
					 ../bin/ssradar \
					 ../bin/sza2time \
					 ../bin/time2sza \
					 ../bin/uvspec\
					 ../bin/uvspecfunction \
					 ../bin/uvspec_mpi \
					 ../bin/vac2air \
					 ../bin/vrs_ocean_albedo \
					 ../bin/zenith \
					 $(CLDGEN_TOOL) \
					 $(WORLDLOOP_TARGETS) \

UNINSTALL_BINARIES = $(patsubst ../bin/, , $(BINARIES))

all: dep $(BINARIES)
.PHONY: all

dep: $(SRC_C) sbtaugas.h
	$(CC) $(CPPFLAGS) $(ALL_CFLAGS) -MM $(SRC_C)  > .depend
.PHONY: dep

../bin/add_worldloop_result: $(ADD_WORLDLOOP_RESULT_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/angres: $(ANGRES_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

cldgen_lex.c: cldgen_lex.l
	$(LEX) $(LFLAGS) -o $@ $<

../bin/cldgen: $(CLDGEN_OBJECTS) cldgen_lex.o $(RADTRAN_C_LIB)
	$(lcommand)
	$(CC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/conv: $(CONVOLUTION_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/integrate: $(INTEGRATE_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/make_slitfunction: $(MAKE_SLITFUNCTION_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/make_angresfunc: $(MAKE_ANGRESFUNC_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

mie_lex.c: mie_lex.l
	$(LEX) $(LFLAGS) -o $@ $<

../bin/mie: $(MIE_OBJECTS) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/noon: $(NOON_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/phase: $(PHASE_OBJECTS) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/pmom: $(PMOM_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/read_Stamnes_tab: $(DER_STAM_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
  $(ALL_LDFLAGS)

../bin/snowalbedo: $(SNOWALBEDO_OBJECTS) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/sofi: $(SOFI_OBJECTS) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/spline: $(SPLINE_INTERPOLATION_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/sza2time: $(SZA2TIME_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/time2sza: $(TIME2SZA_OBJECTS) $(RADTRAN_F_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

uvspec_lex.c: $(PYTHONDIR)uvspec_lex.l
	$(LEX) $(LFLAGS) -o $@ $<


../bin/uvspec: $(UVSPEC_OBJECTS) uvspec_lex.o $(RADTRAN_LIB)
	make -e -C ../src_py || exit $$?; cd ../src # Not very elegant, but it seems to do stuff when *py files are updated aky20120804
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/arlem: $(ARLEM_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/uvspecfunction: $(UVSPEC_OBJECTS) $(UVSPEC_FUNCTION_OBJECTS) uvspecfunction_lex.o $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/uvspec_mpi: $(UVSPEC_OBJECTS) $(UVSPEC_MPI_OBJECTS) uvspecfunction_lex.o $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/vrs_ocean_albedo: $(VRSALBEDO_OBJECTS) $(UVSPEC_OBJECTS) $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/worldloop: $(UVSPEC_OBJECTS) worldloop.o uvspecfunction_lex.o $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/check_worldloop_result: check_worldloop_result.o $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/create_worldloop_grid: create_worldloop_grid.o $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/copy_worldloop_grid: copy_worldloop_grid.o $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/cldprp: $(CLDPRP_OBJECTS) $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/plkavg: $(PLKAVG_OBJECTS) $(RADTRAN_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $+ \
	$(ALL_LDFLAGS)

../bin/zenith: $(ZENITH_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ $^ \
	$(ALL_LDFLAGS)

../bin/vac2air: $(VAC2AIR_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ \
	$(VAC2AIR_OBJECTS) \
	$(RADTRAN_C_LIB) \
	$(ALL_LDFLAGS)

../bin/ssradar: $(SSRADAR_OBJECTS) $(RADTRAN_C_LIB)
	$(lcommand)
	$(CCC) $(CFLAGS) -o $@ \
	$(SSRADAR_OBJECTS) \
	$(RADTRAN_C_LIB) $(ALL_LDFLAGS)

sbtaugas.inc: sbtaugas.fparam sbtaugas.param
	cp sbtaugas.fparam $@
	$(AWK) 'NF==2{print "      parameter ("$$1" = "$$2")"}' sbtaugas.param >> $@

sbtaugas.h: sbtaugas.param
	$(AWK) 'NF==2{print "#define "$$1" "$$2}' $< > $@

## -lccmalloc -ldl

# How to make .o files from .c files. We need to modify the rule because
# we don't have header files here in contrast to the src_c directory.
# We also re-compile if the libraries have changed, assuming that those have
# been recompiled if the header files in src_c and src_f have been changed;
# this is not very nice and we should rather check which source file actually
# includes a header from src_c/src_f and only recompile those
%.o: %.c $(UVSPEC_HEADER) Makefile ../Makeconf
	$(ccommand)
	$(CC) -c $(CPPFLAGS) $(ALL_CFLAGS) $<

# special rule for uvspecfunction_lex.o because we need to #define UVSPEC_FUNCTION
uvspecfunction_lex.o: uvspec_lex.c $(UVSPEC_HEADER) $(RADTRAN_F_LIB) $(RADTRAN_C_LIB)
	$(ccommand)
	$(CC) -o $@ -c $(CPPFLAGS) $(ALL_CFLAGS) -DUVSPEC_FUNCTION=1 $<

# special dependencies:
ck.o       : $(LOWTRAN_CINCLUDE)
taugas.o   : $(LOWTRAN_FINCLUDE) Makefile ../Makeconf
sbtaugas.o : $(LOWTRAN_FINCLUDE) Makefile ../Makeconf
taucor.o   : $(LOWTRAN_FINCLUDE) Makefile ../Makeconf
depthscl.o : $(LOWTRAN_FINCLUDE) Makefile ../Makeconf

$(patsubst %.f,%.o,$(FULIOU_FSOURCES)) : Makefile ../Makeconf $(FULIOU_INCLUDES)
$(patsubst %.c,%.o,$(FULIOU_CSOURCES)) : Makefile ../Makeconf $(FULIOU_INCLUDES)

check: all
.PHONY: check

install:
	$(MKPATH) $(bindir)
	$(INSTALL) $(BINARIES) $(bindir)
.PHONY: install

uninstall:
	for f in $(BINARIES); do \
	  fbase=`basename $$f`; \
	  $(RM) $(bindir)/$$fbase; \
	done
.PHONY: uninstall

clean:
	$(RM) -f *.a *.o .depend *.gcno *.gcda
.PHONY: clean

mostlyclean:
	$(RM) -f *.o
.PHONY: mostlyclean

distclean: clean
	$(RM) -f Makefile
.PHONY: distclean

realclean: distclean
	$(RM) -f Makefile
	$(RM) $(BINARIES)
.PHONY: realclean

dist:
	touch lex.c
	ln $(DISTFILES) ../`cat ../.fname`/src
.PHONY: dist

mystic-dist:
	touch lex.c
	ln $(MYSTICDISTFILES) ../`cat ../.fname`/src
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status

-include .depend
