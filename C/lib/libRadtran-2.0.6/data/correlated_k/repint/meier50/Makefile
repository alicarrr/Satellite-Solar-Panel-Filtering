#
# Makefile for libRadtran's data/correlated_k/repint/meier50 directory
#


TOPDIR = ../../../../

srcdir = .


THISDIR = meier50

DATAFILES =  $(wildcard *.cdf) 

DISTFILES = Makefile.in $(DATAFILES)

MYSTICDISTFILES = $(DISTFILES)

include $(TOPDIR)/Makeconf

all: 

dist: all
	ln $(DISTFILES) ../../../../`cat ../../../../.fname`/data/correlated_k/repint/$(THISDIR)
.PHONY: dist

mystic-dist: all
	ln $(MYSTICDISTFILES) ../../../../`cat ../../../../.fname`/data/correlated_k/repint/$(THISDIR)
.PHONY: mystic-dist

install: all
	$(MKPATH) $(datafiledir)/correlated_k/repint/$(THISDIR)
	$(INSTALL_DATA) $(DATAFILES) $(datafiledir)/correlated_k/repint/$(THISDIR)
.PHONY: install

clean:

mostlyclean: clean

distclean: clean
	rm -f Makefile config.status

realclean: distclean

Makefile: Makefile.in ../../../../config.status
	cd ../../../../; sh config.status

# Prevent GNU make v3 from overflowing arg limit on SysV.
.NOEXPORT:

# eof
