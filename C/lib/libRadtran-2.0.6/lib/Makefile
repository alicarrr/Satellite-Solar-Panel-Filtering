#
# Makefile for libRadtran's lib directory
#
# Arve Kylling
# arve.kylling@nilu.no

TOPDIR = ..

include $(TOPDIR)/Makeconf

srcdir = .

WHOLE_ARCHIVE = -Wl,--whole-archive
NO_WHOLE_ARCHIVE = -Wl,--no-whole-archive

DISTFILES = Makefile.in

MYSTICDISTFILES = $(DISTFILES)

all: libRadtran.so 
.PHONY: all

libRadtran.so: libRadtran_c.a libRadtran_f.a
	$(CC) $(CFLAGS) --shared -o $@ $(WHOLE_ARCHIVE) libRadtran_c.a $(NO_WHOLE_ARCHIVE)  libRadtran_f.a $(FLIBS) $(LEXLIB) $(LIBNETCDF) $(LIBGSL) $(LIBGMP) $(LIBOPENGL) $(ALL_LDFLAGS) -lm

install: all
	$(MKPATH) $(libdir)
	$(INSTALL_DATA) libRadtran.so $(libdir)/libRadtran.so
	$(INSTALL_DATA) LibRadtranConfig.cmake.install $(libdir)/LibRadtranConfig.cmake
.PHONY: install

uninstall:
.PHONY: uninstall

spell: $(SPELL)
.PHONY: spell

clean:
	rm libRadtran.so
.PHONY: clean

mostlyclean:
.PHONY: mostlyclean

distclean: clean
	rm -f Makefile
.PHONY: distclean

realclean: distclean
	rm -f tags TAGS
.PHONY: realclean

dist:
	ln -f $(DISTFILES) ../`cat ../.fname`/lib
.PHONY: dist

mystic-dist:
	ln -f $(MYSTICDISTFILES) ../`cat ../.fname`/lib
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status

