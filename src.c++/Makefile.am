## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rheinfall.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/mathdefs.hpp

# build example programs
EXTRA_PROGRAMS = irank-omp irank-mpi irank-mpi-omp
bin_PROGRAMS = irank
if HAVE_OPENMP_CXX
  bin_PROGRAMS +=  irank-omp
endif
if HAVE_MPICXX
  bin_PROGRAMS += irank-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += irank-mpi-omp
endif
endif

irank_SOURCES = examples/irank.cpp
irank_omp_SOURCES = examples/irank.cpp
irank_mpi_SOURCES = examples/irank.cpp
irank_mpi_omp_SOURCES = examples/irank.cpp

AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c++/rheinfall $(BOOST_CPPFLAGS)
AM_CXXFLAGS = 

irank_omp_CPPFLAGS= $(AM_CPPFLAGS) $(OPENMP_CPPFLAGS)
irank_omp_CXXFLAGS= $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
irank_omp_LDFLAGS= # openmp?
irank_omp_LIBS= # openmp?

irank_mpi_CPPFLAGS= $(AM_CPPFLAGS) -DMPI_SERIALIZED
irank_mpi_CXXFLAGS= $(AM_CXXFLAGS)
irank_mpi_LDFLAGS= $(BOOST_SERIALIZATION_LDFLAGS)
irank_mpi_LIBS= $(BOOST_SERIALIZATION_LIBS)

irank_mpi_omp_CPPFLAGS= $(irank_omp_CPPFLAGS) $(irank_mpi_CPPFLAGS)
irank_mpi_omp_CXXFLAGS= $(irank_omp_CXXFLAGS) $(irank_mpi_CXXFLAGS)
irank_mpi_omp_LDFLAGS= $(irank_omp_LDFLAGS) $(irank_mpi_LDFLAGS)
irank_mpi_omp_LIBS= $(irank_omp_LIBS) $(irank_mpi_LIBS)


# unit tests
check_PROGRAMS = test-sparserow test-denserow

test_sparserow_SOURCES = test/test_sparserow.cpp
test_sparserow_CXXFLAGS = $(BOOST_CPPFLAGS)
test_sparserow_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_sparserow_LIBS= $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

test_denserow_SOURCES = test/test_denserow.cpp
test_denserow_CXXFLAGS = $(BOOST_CPPFLAGS)
test_denserow_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_denserow_LIBS= $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS = $(check_PROGRAMS)
