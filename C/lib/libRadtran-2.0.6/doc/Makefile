#
# Makefile for libRadtran's doc directory
#
# Arve Kylling
# arve.kylling@itek.norut.no

TOPDIR = ..

include $(TOPDIR)/Makeconf

srcdir = .
VPATH = ../src/

SUBDIRS = figs

F_SOURCES = ../libsrc_f/wcloud.f
F_DIR     = ../libsrc_f/

# Currently not used
#C_SOURCES = ../libsrc_c/ascii.c ../libsrc_c/cnv.c ../libsrc_c/equation.c \
#            ../libsrc_c/function.c ../libsrc_c/integrat.c ../libsrc_c/linear.c \
#            ../libsrc_c/regress.c ../libsrc_c/spl.c ../libsrc_c/sun.c \
#	    ../libsrc_c/fortran_and_c.c ../libsrc_c/wavelength.c \
#            ../libsrc_c/snow.c ../libsrc_c/miecalc.c $(CLDGEN_CSOURCES)

C_DIR     = ../libsrc_c/

doc: libRadtran.pdf

#------------------------------------------------------
# Documentation of options (extracted from .lex files)
#------------------------------------------------------

L_SOURCES = ../src/mie_lex.l ../src_py/uvspec_lex.l
L_LATEX   = $(patsubst ../src_py/%, %, $(patsubst ../src/%, %, $(patsubst %.l, %.tex, $(L_SOURCES))))
L_PAIRS  := $(join ${L_LATEX},$(addprefix :,${L_SOURCES}))
$(foreach r,$(L_PAIRS),$(eval $r))

$(L_LATEX):
	gawk -f sdoc2.awk < $< > $@

#------------------------------------------------------
# Documentation of tools (extracted from .c files)
#------------------------------------------------------

# Please include source code file from which documentation should go into
# the user guide here
T_SOURCES = \
						../src/angres.c \
						../src/make_angresfunc.c \
						../src/make_slitfunction.c \
						../src/noon.c	\
						../src/phase.c \
						../src/pmom.c \
						../src/read_Stamnes_tab.c \
						../src/snowalbedo.c \
						../src/zenith.c	\

# Commented for ESA distribution ../src/worldloop.c

T_LATEX  = $(patsubst ../src/%, %, $(patsubst %.c, %.tex, $(T_SOURCES)))
T_PAIRS := $(join ${T_LATEX},$(addprefix :,${T_SOURCES}))
$(foreach r,$(T_PAIRS),$(eval $r))

$(T_LATEX):
	gawk -f sdoc2.awk < $< > $@


#-----------------------------------------------------------
# Other files (not automatically extracted from source code)
#-----------------------------------------------------------

LATEX  = \
	intro.tex \
	libRadtran.tex \
	mystic.tex \
	optical_properties.tex \
	radiative_transfer.tex \
	radiative_transfer_theory.tex \
	ssradar.tex \
	title.tex \
	tools.tex \

STY = cprotect.dtx cprotect.sty

BIBTEX = literature.bib egu.bst

FORMATTED = libRadtran.pdf

IMAGES =

#20100312 All images are included by Makefile in figs directory
#./figs/figure3.pdf ./figs/fuyang.pdf ./figs/lblo2a.pdf ./figs/sbdarto2a.pdf\
#	 ./figs/sbdart.pdf ./figs/transpz.pdf ./figs/mystic2.png


#----------------------------------------------------------
# Distribution files
#----------------------------------------------------------

DISTFILES = Makefile.in $(LATEX) $(STY) $(BIBTEX) $(FORMATTED) $(IMAGES) sdoc2.awk  # README.MC

MYSTICDISTFILES = $(DISTFILES)

#---------------------------------------------------------
# Target specifications
#---------------------------------------------------------

all: libRadtran.pdf
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: all

undocumented:
	echo Error, currently not implemented!
#	cp mystic.sdoc mystic.sdoc.save
#	echo @set mystic >> mystic.sdoc
#	echo @set undocumented >> mystic.sdoc
#	$(MAKE) all
#	mv mystic.sdoc.save mystic.sdoc

.PHONY: all

libRadtran.pdf: $(LATEX) $(BIBTEX) $(L_LATEX) $(T_LATEX)
	echo \\setboolean{mystic_doc}{true} > mystic.doc
	echo \\setboolean{mystic3d_doc}{false} >> mystic.doc
	echo \\setboolean{lidar_doc}{false} >> mystic.doc
	pdflatex libRadtran.tex
	bibtex libRadtran; pdflatex libRadtran.tex; pdflatex libRadtran.tex;
	makeindex libRadtran
	pdflatex libRadtran.tex; pdflatex libRadtran.tex

#AK 20100308 We do not want install to depend on all as then the documentation is
# rebuilt everytime a make install is performed.
install:
	$(MKPATH) $(docdir)
	$(INSTALL) -m644 libRadtran.pdf $(docdir)
.PHONY: install

uninstall:
	$(RM) $(docdir)/libRadtran.pdf
.PHONY: uninstall

spell: $(SPELL)
.PHONY: spell

clean:
	rm -f libRadtran.aux libRadtran.log libRadtran.toc \
	libRadtran.bbl libRadtran.blg libRadtran.idx libRadtran.ilg libRadtran.ind libRadtran.out \
  tools.aux radiative_transfer.aux optical_properties.aux radiative_transfer_theory.aux \
	uvspec_lex.aux intro.aux mie_lex.aux \
	$(L_LATEX) $(T_LATEX)

.PHONY: clean

mostlyclean:
.PHONY: mostlyclean

distclean: clean
	rm -f Makefile
.PHONY: distclean

realclean: distclean
	rm -f tags TAGS
.PHONY: realclean

dist: all
	for dir in $(SUBDIRS); do mkdir ../`cat ../.fname`/doc/$$dir; $(MAKE) -C $$dir $@ || exit $$?; done
	ln $(DISTFILES) ../`cat ../.fname`/doc
.PHONY: dist

mystic-dist: all
	for dir in $(SUBDIRS); do mkdir ../`cat ../.fname`/doc/$$dir; $(MAKE) -C $$dir $@ || exit $$?; done
	ln $(MYSTICDISTFILES) ../`cat ../.fname`/doc
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status

list-lexfiles:
	$(info L_SOURCES:$(L_SOURCES))
	$(info L_LATEX  :$(L_LATEX))
	$(info T_SOURCES:$(T_SOURCES))
	$(info T_LATEX  :$(T_LATEX))
.PHONY: list-lexfiles
