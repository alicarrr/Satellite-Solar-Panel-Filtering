TOPDIR = ..
include $(TOPDIR)/Makeconf

DISTFILES = \
	catch.hpp \
	main.cpp \
	Makefile.in \
	test_allocnd.cpp \
	test_common_math.cpp \
	test_triangle_surface.cpp \
	test_tribox3.cpp

MYSTICDISTFILES = $(DISTFILES)

RADTRAN_F_LIB =  $(TOPDIR)/lib/libRadtran_f.a
RADTRAN_C_LIB =  $(TOPDIR)/lib/libRadtran_c.a

# duplicate for the linker to resolve circular dependencies
RADTRAN_LIB = $(RADTRAN_C_LIB) $(RADTRAN_F_LIB) $(RADTRAN_C_LIB) $(RADTRAN_F_LIB)

objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

CXXFLAGS := ${CXXFLAGS} -I..

all: exec_tests
.PHONY: all

run_unittests: exec_tests
	./exec_tests
.PHONY: run_unittests

exec_tests: ${objects} ${RADTRAN_LIB}
	${CXX} ${ALL_CFLAGS} ${CXXFLAGS} -o $@ $+ \
	$(ALL_LDFLAGS)

%.o: %.cpp
	${CXX} ${ALL_CFLAGS} ${CXXFLAGS} -c $<

clean:
	$(RM) -f ${objects} exec_tests
.PHONY: clean

dist:
	ln $(DISTFILES) ../`cat ../.fname`/unittests
.PHONY: dist

mystic-dist:
	ln $(MYSTICDISTFILES) ../`cat ../.fname`/unittests
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status
