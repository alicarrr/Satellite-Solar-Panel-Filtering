#
# Makefile for Radtran
#
# Arve Kylling
# arve.kylling@gmail.com
#

TOPDIR = .

srcdir = .


include Makeconf

DISTFILES = COPYING ChangeLog INSTALL README README.unpacking\
        README.windows README.macosx \
        TODO\
	Makefile.in Makeconf.in\
	configure configure.in config.guess config.sub install-sh

MYSTICDISTFILES = $(DISTFILES)

# Complete directory trees to distribute.
DISTDIRS = m4

# Subdirectories to run make in for the primary targets.
SUBDIRS = examples src_py libsrc_c libsrc_f src lib unittests GUI

# Subdirectories to run make in for the installation.
INSTALLSUBDIRS = src_py libsrc_c libsrc_f src examples doc data GUI

# Subdirectories to run `make dist' in
DISTSUBDIRS = $(SUBDIRS) data sdoc doc flexstor test bin

# Subdirectories to run `make cdisort-dist' in
CDISORTDISTSUBDIRS = libsrc_c

all:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: all

uvspec_lex:
	$(MAKE) -e -C src_py uvspec_lex.l
.PHONY: uvspec_lex

libRadtran_c:
	$(MAKE) -e -C libsrc_c
.PHONY: libRadtran_c

libRadtran_f:
	$(MAKE) -e -C libsrc_f
.PHONY: libRadtran_f

arlem: libRadtran_c libRadtran_f
	$(MAKE) -e -C src ../bin/arlem
.PHONY: arlem

uvspec: uvspec_lex libRadtran_c libRadtran_f
	$(MAKE) -e -C src ../bin/uvspec
.PHONY: uvspec

src: uvspec_lex libRadtran_c libRadtran_f
	$(MAKE) -e -C src
.PHONY: src

check: all
	$(MAKE) -e -C test
.PHONY: check

unittests: libRadtran_c libRadtran_f
	$(MAKE) -e -C unittests
.PHONY: unittests

run_unittests: unittests
	$(MAKE) -e -C unittests run_unittests
.PHONY: run_unittests

verbosecheck: all
	$(MAKE) -e -C test verbosecheck
.PHONY: verbosecheck

valgrind: all
	$(MAKE) -e -C test valgrind
.PHONY: verbosecheck

globalcheck: all
	$(MAKE) -e -C test globalcheck
.PHONY: globalcheck

answers: all
	$(MAKE) -e -C test answers
.PHONY: answers

checkeverything: all
	$(MAKE) -e -C test everything
.PHONY: checkeverything

doc: uvspec_lex
	$(MAKE) -e -C doc
.PHONY: doc

install:
	for dir in $(INSTALLSUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: install

uninstall:
	for dir in $(INSTALLSUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
	$(RM) -rf $(datadir)

#	$(RM) -rf $(datadir)/examples; \
#	$(RM) -rf $(datadir)/doc;

.PHONY: uninstall

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: clean

mostlyclean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: mostlyclean

distclean: clean
	rm -f *.a Makefile Makeconf config.status
.PHONY: distclean

realclean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
	rm -f *.a Makefile Makeconf config.status
.PHONY: realclean

links-for-dist:
	echo libRadtran-$(version) > .fname
	rm -rf `cat .fname`
	mkdir `cat .fname`
	ln $(DISTFILES) `cat .fname`
	for dir in $(DISTDIRS); do cp -al $$dir `cat .fname`/; done
	for dir in $(DISTSUBDIRS); do \
			mkdir `cat .fname`/$$dir; \
			cd $$dir; \
			$(MAKE) dist || exit $$?; \
			cd ..; \
		done
.PHONY: links-for-dist

dist: links-for-dist
	tar chf `cat .fname`.tar `cat .fname`
	rm -rf `cat .fname`
	tar xf `cat .fname`.tar
	chmod -R a+rwX `cat .fname`
	tar cf `cat .fname`.tar `cat .fname`
	rm -rf `cat .fname` .fname
.PHONY: dist

links-for-mystic-dist:
	echo libRadtran-$(mystic-version) > .fname
	rm -rf `cat .fname`
	mkdir `cat .fname`
	ln $(MYSTICDISTFILES) `cat .fname`
	for dir in $(DISTSUBDIRS); do mkdir `cat .fname`/$$dir; cd $$dir; $(MAKE) mystic-dist || exit $$?; cd ..; done
.PHONY: links-for-mystic-dist

mystic-dist: links-for-mystic-dist
	tar chf `cat .fname`.tar `cat .fname`
	rm -rf `cat .fname`
	tar xf `cat .fname`.tar
	chmod -R a+rwX `cat .fname`
	tar cf `cat .fname`.tar `cat .fname`
	rm -rf `cat .fname` .fname
.PHONY: mystic-dist

dist-z:
	echo Radtran-$(version).tar > .fname
	gzip --best `cat .fname`
.PHONY: dist-z

# Make distribution for cdisort etc.

links-for-cdisort-dist:
	echo cdisort-$(cdisort-version) > .cname
	rm -rf `cat .cname`
	mkdir `cat .cname`
	for dir in $(CDISORTDISTSUBDIRS); do cd $$dir; $(MAKE) cdisort-dist || exit $$?; cd ..; done
.PHONY: links-for-cdisort-dist

cdisort-dist: links-for-cdisort-dist
	tar chf `cat .cname`.tar `cat .cname`
	rm -rf `cat .cname`
	tar xf `cat .cname`.tar
	chmod -R a+rwX `cat .cname`
	tar cf `cat .cname`.tar `cat .cname`
	rm -rf `cat .cname` .cname
.PHONY: cdisort-dist

Makefile: Makefile.in ./config.status
	sh config.status
