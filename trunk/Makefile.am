## Process this file with automake to produce Makefile.in
bin_PROGRAMS = irank
irank_SOURCES = src.c/main.c src.c/rheinfall.c src.c/vpu.c src.c/switchboard.c src.c/comm.c src.c/row.c src.c/ge.c 


# FIXME: GCC-specific flags!
AM_CFLAGS = -Wall -DWITH_MPI
SUBDIRS =
LDADD =


# GNUlib support.  Output from `gnulib-tool --import xalloc xalloc-die`:
#  - mention "lib" in SUBDIRS in Makefile.am,
#  - mention "-I m4" in ACLOCAL_AMFLAGS in Makefile.am,
#  - mention "m4/gnulib-cache.m4" in EXTRA_DIST in Makefile.am,
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = m4/gnulib-cache.m4
AM_CFLAGS += -I$(top_builddir)/lib -I$(top_srcdir)/lib
LDADD += lib/libgnu.a
SUBDIRS += lib


# tests
TESTS = test-xarray
check_PROGRAMS = test-xarray
test_xarray_SOURCES = src.c/test-xarray.c
