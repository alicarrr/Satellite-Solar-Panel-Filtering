#
# Makefile for generating uvspec_lex.l from src_py dictionary
#

TOPDIR = ..

include $(TOPDIR)/Makeconf

PYTHONFILES = \
	GUI_definition.py \
	aerosol_options.py \
	cloud_options.py \
	general_atmosphere_options.py \
	geometry_options.py \
	lex_starter.l \
	mc_options.py \
	molecular_options.py \
	option_definition.py \
	output_options.py \
	profile_options.py \
	solver_options.py \
	special_options.py \
	spectral_options.py \
	surface_options.py \
	translate.py \
	writeLex.py.in \

DISTFILES = Makefile.in $(PYTHONFILES) uvspec_lex.l writeLex.py

MYSTICDISTFILES = $(DISTFILES)

all: uvspec_lex.l
.PHONY: all

uvspec_lex.l: $(PYTHONFILES) writeLex.py
	$(PYTHON) writeLex.py

clean:
.PHONY: clean

distclean: clean
	rm -f Makefile
.PHONY: distclean

realclean: distclean
.PHONY: realclean

install:
.PHONY: install

uninstall:
.PHONY: uninstall

dist:
	ln $(DISTFILES) ../`cat ../.fname`/src_py
.PHONY: dist

mystic-dist:
	ln $(MYSTICDISTFILES) ../`cat ../.fname`/src_py
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status

writeLex.py: writeLex.py.in
	cd ..; sh config.status
