#
# Makefile for libRadtran's example/<subdir> directory
#

TOPDIR = ../..

srcdir = .


THISDIR = $(notdir $(CURDIR))

SUBDIRS =

EXAMPLEFILES = \
	$(wildcard *.dat) \
	$(wildcard *.DAT) \
	$(wildcard *.INP) \
	$(wildcard *.OUT) \

DISTFILES = Makefile.in $(EXAMPLEFILES)

include $(TOPDIR)/Makeconf

all:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ || exit $$?; done
.PHONY: all

dist: all
	for dir in $(SUBDIRS); do mkdir $(TOPDIR)/`cat $(TOPDIR)/.fname`/examples/$(THISDIR)/$$dir; $(MAKE) -C $$dir $@ || exit $$?; done; \
	ln $(DISTFILES) $(TOPDIR)/`cat $(TOPDIR)/.fname`/examples/$(THISDIR)
.PHONY: dist

mystic-dist: dist
.PHONY: mystic-dist

install: all
	$(MKPATH) $(examplesdir)/$(THISDIR)
	$(INSTALL_DATA) $(EXAMPLEFILES) $(examplesdir)/$(THISDIR)
.PHONY: install

clean:
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
