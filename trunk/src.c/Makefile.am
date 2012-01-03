## Process this file with automake to produce Makefile.in

SUBDIRS = gnulib

# binaries built here
bin_PROGRAMS = crank-int crank-double
if HAVE_MPICC
  bin_PROGRAMS += crank-int-mpi crank-double-mpi 
endif

AM_CPPFLAGS =
if WITH_GE
  AM_CPPFLAGS += -DWITH_GE
endif
AM_CFLAGS =
AM_LDFLAGS =
LDADD =

crank_int_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
crank_int_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c/gnulib -DWITH_INT_VALUES
crank_int_CFLAGS = $(AM_CFLAGS)
crank_int_LDFLAGS = $(AM_LDFLAGS)
crank_int_LDADD = $(LDADD)

crank_int_mpi_SOURCES = $(crank_int_SOURCES)
crank_int_mpi_CPPFLAGS = $(crank_int_CPPFLAGS) -DWITH_MPI $(MPI_CFLAGS)
crank_int_mpi_CFLAGS = $(crank_int_CFLAGS) 
crank_int_mpi_LDFLAGS = $(crank_int_CFLAGS) $(MPI_CLDFLAGS)
crank_int_mpi_LDADD = $(crank_int_LDADD) 


crank_double_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
crank_double_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c/gnulib -DWITH_DOUBLE_VALUES
crank_double_CFLAGS = $(AM_CFLAGS)
crank_double_LDFLAGS = $(AM_LDFLAGS)
crank_double_LDADD = $(LDADD)

crank_double_mpi_SOURCES = $(crank_double_SOURCES)
crank_double_mpi_CPPFLAGS = $(crank_double_CPPFLAGS) -DWITH_MPI $(MPI_CFLAGS)
crank_double_mpi_CFLAGS = $(crank_double_CFLAGS) 
crank_double_mpi_LDFLAGS = $(crank_double_CFLAGS) $(MPI_CLDFLAGS)
crank_double_mpi_LDADD = $(crank_double_LDADD) 


# GNUlib support.  Output from `gnulib-tool --import xalloc xalloc-die`:
#  - mention "lib" in SUBDIRS in Makefile.am,
#  - mention "-I m4" in ACLOCAL_AMFLAGS in Makefile.am,
#  - mention "m4/gnulib-cache.m4" in EXTRA_DIST in Makefile.am,
ACLOCAL_AMFLAGS = -I$(top_srcdir)/build-aux/m4
EXTRA_DIST = $(top_srcdir)/build_aux/m4/gnulib-cache.m4
AM_CFLAGS += -I$(builddir)/gnulib -I$(srcdir)/gnulib
LDADD += $(builddir)/gnulib/libgnu.a


# tests
check_PROGRAMS = test-xarray
test_xarray_SOURCES = test-xarray.c

TESTS = test-xarray
