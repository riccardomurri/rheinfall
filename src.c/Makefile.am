## Process this file with automake to produce Makefile.in

SUBDIRS = gnulib

# binaries built here
EXTRA_PROGRAMS = irank-homo-omp irank-homo-mpi irank-homo-mpi-omp
bin_PROGRAMS = irank-homo
if HAVE_OPENMP_CC
  bin_PROGRAMS +=  irank-homo-omp
endif
if HAVE_MPICC
  bin_PROGRAMS += irank-homo-mpi
if HAVE_OPENMP_CC
    bin_PROGRAMS += irank-homo-mpi-omp
endif
endif

irank_homo_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
AM_CFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c/gnulib
LDADD =

irank_homo_omp_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
irank_homo_omp_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS)
irank_homo_omp_LDFLAGS = # openmp?
irank_homo_omp_LIBS = # openmp?

irank_homo_mpi_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
irank_homo_mpi_CFLAGS = $(AM_CFLAGS) -DMPI_SERIALIZED
irank_homo_mpi_LDFLAGS =
irank_homo_mpi_LIBS =

irank_homo_mpi_omp_SOURCES = main.c rheinfall.c vpu.c switchboard.c comm.c row.c ge.c 
irank_homo_mpi_omp_CFLAGS = $(irank_homo_omp_CFLAGS) $(irank_homo_mpi_CFLAGS)
irank_homo_mpi_omp_LDFLAGS = $(irank_homo_omp_LDFLAGS) $(irank_homo_mpi_LDFLAGS)
irank_homo_mpi_omp_LIBS = $(irank_homo_omp_LIBS) $(irank_homo_mpi_LIBS)


# GNUlib support.  Output from `gnulib-tool --import xalloc xalloc-die`:
#  - mention "lib" in SUBDIRS in Makefile.am,
#  - mention "-I m4" in ACLOCAL_AMFLAGS in Makefile.am,
#  - mention "m4/gnulib-cache.m4" in EXTRA_DIST in Makefile.am,
ACLOCAL_AMFLAGS = -I$(top_srcdir)/build-aux/m4
EXTRA_DIST = $(top_srcdir)/build_aux/m4/gnulib-cache.m4
AM_CFLAGS += -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
LDADD += gnulib/libgnu.a


# tests
check_PROGRAMS = test-xarray
test_xarray_SOURCES = test-xarray.c

TESTS = test-xarray
