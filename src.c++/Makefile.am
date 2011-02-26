## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rheinfall.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/mathdefs.hpp

# build example programs
bin_PROGRAMS = irank drank qrank zrank
EXTRA_PROGRAMS = irank-omp drank-omp qrank-omp zrank-omp
if HAVE_OPENMP_CXX
  bin_PROGRAMS +=  irank-omp drank-omp qrank-omp zrank-omp
endif
EXTRA_PROGRAMS += irank-mpi drank-mpi qrank-mpi zrank-mpi
EXTRA_PROGRAMS += irank-mpi-omp drank-mpi-omp qrank-mpi-omp zrank-mpi-omp
if HAVE_MPICXX
  bin_PROGRAMS += irank-mpi drank-mpi #qrank-mpi zrank-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += irank-mpi-omp drank-mpi-omp #qrank-mpi-omp zrank-mpi-omp
endif
endif

if HAVE_OPENMP_CXX
omp: irank-omp drank-omp qrank-omp zrank-omp
endif
if HAVE_MPICXX
mpi: irank-mpi drank-mpi #qrank-mpi zrank-mpi
if HAVE_OPENMP_CXX
mpi-omp: irank-mpi-omp drank-mpi-omp #qrank-mpi-omp zrank-mpi-omp
endif
endif


AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c++/rheinfall $(BOOST_CPPFLAGS)
AM_CXXFLAGS = 

# irank - example program using `long long` integer values
irank_SOURCES = examples/irank.cpp
irank_CPPFLAGS= $(AM_CPPFLAGS) 
irank_CXXFLAGS= $(AM_CXXFLAGS)
irank_LDFLAGS= $(AM_LDFLAGS) 
irank_LDADD = $(AM_LIBS) 

irank_omp_SOURCES = $(irank_SOURCES)
irank_omp_CPPFLAGS = $(irank_CPPFLAGS) $(OPENMP_CPPFLAGS)
irank_omp_CXXFLAGS = $(irank_CXXFLAGS) $(OPENMP_CXXFLAGS)
irank_omp_LDFLAGS = $(irank_LDFLAGS) # openmp?
irank_omp_LDADD = $(irank_LDADD) # openmp?

irank_mpi_SOURCES = $(irank_SOURCES)
irank_mpi_CPPFLAGS = $(irank_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
irank_mpi_CXXFLAGS = $(irank_CXXFLAGS) 
irank_mpi_LDFLAGS = $(irank_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
irank_mpi_LDADD = $(irank_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

irank_mpi_omp_SOURCES = $(irank_mpi_SOURCES)
irank_mpi_omp_CPPFLAGS= $(irank_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
irank_mpi_omp_CXXFLAGS= $(irank_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
irank_mpi_omp_LDFLAGS= $(irank_mpi_LDFLAGS) # openmp?
irank_mpi_omp_LDADD = $(irank_mpi_LDADD) # openmp?

# drank - example program using `double` floating-point values
drank_SOURCES = examples/drank.cpp
drank_CPPFLAGS= $(AM_CPPFLAGS) 
drank_CXXFLAGS= $(AM_CXXFLAGS)
drank_LDFLAGS= $(AM_LDFLAGS) 
drank_LDADD = $(AM_LIBS) 

drank_omp_SOURCES = $(drank_SOURCES)
drank_omp_CPPFLAGS = $(drank_CPPFLAGS) $(OPENMP_CPPFLAGS)
drank_omp_CXXFLAGS = $(drank_CXXFLAGS) $(OPENMP_CXXFLAGS)
drank_omp_LDFLAGS = $(drank_LDFLAGS) # openmp?
drank_omp_LDADD = $(drank_LDADD) # openmp?

drank_mpi_SOURCES = $(drank_SOURCES)
drank_mpi_CPPFLAGS = $(drank_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
drank_mpi_CXXFLAGS = $(drank_CXXFLAGS)
drank_mpi_LDFLAGS = $(drank_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
drank_mpi_LDADD = $(drank_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

drank_mpi_omp_SOURCES = $(drank_mpi_SOURCES)
drank_mpi_omp_CPPFLAGS= $(drank_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
drank_mpi_omp_CXXFLAGS= $(drank_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
drank_mpi_omp_LDFLAGS= $(drank_mpi_LDFLAGS) # openmp?
drank_mpi_omp_LDADD = $(drank_mpi_LDADD) # openmp?

# qrank - example program using GMP `mpq` infinite precision rational values
qrank_SOURCES = examples/qrank.cpp
qrank_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS)
qrank_CXXFLAGS= $(AM_CXXFLAGS)
qrank_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
qrank_LDADD= $(AM_LIBS) $(GMPXX_LIBS)

qrank_omp_SOURCES = $(qrank_SOURCES)
qrank_omp_CPPFLAGS = $(qrank_CPPFLAGS) $(OPENMP_CPPFLAGS)
qrank_omp_CXXFLAGS = $(qrank_CXXFLAGS) $(OPENMP_CXXFLAGS)
qrank_omp_LDFLAGS = $(qrank_LDFLAGS) # openmp?
qrank_omp_LDADD = $(qrank_LDADD) # openmp?

qrank_mpi_SOURCES = $(qrank_SOURCES)
qrank_mpi_CPPFLAGS = $(qrank_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
qrank_mpi_CXXFLAGS = $(qrank_CXXFLAGS)
qrank_mpi_LDFLAGS = $(qrank_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
qrank_mpi_LDADD = $(qrank_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

qrank_mpi_omp_SOURCES = $(qrank_mpi_SOURCES)
qrank_mpi_omp_CPPFLAGS= $(qrank_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
qrank_mpi_omp_CXXFLAGS= $(qrank_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
qrank_mpi_omp_LDFLAGS= $(qrank_mpi_LDFLAGS) # openmp?
qrank_mpi_omp_LDADD = $(qrank_mpi_LDADD) # openmp?

# zrank - example program using GMP `mpz_t` infinite-precision integer values
zrank_SOURCES = examples/zrank.cpp
zrank_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS)
zrank_CXXFLAGS= $(AM_CXXFLAGS)
zrank_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
zrank_LDADD = $(AM_LIBS) $(GMPXX_LIBS)

zrank_omp_SOURCES = $(zrank_SOURCES)
zrank_omp_CPPFLAGS = $(zrank_CPPFLAGS) $(OPENMP_CPPFLAGS)
zrank_omp_CXXFLAGS = $(zrank_CXXFLAGS) $(OPENMP_CXXFLAGS)
zrank_omp_LDFLAGS = $(zrank_LDFLAGS) # openmp?
zrank_omp_LDADD = $(zrank_LDADD) # openmp?

zrank_mpi_SOURCES = $(zrank_SOURCES)
zrank_mpi_CPPFLAGS = $(zrank_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
zrank_mpi_CXXFLAGS = $(zrank_CXXFLAGS)
zrank_mpi_LDFLAGS = $(zrank_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
zrank_mpi_LDADD = $(zrank_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

zrank_mpi_omp_SOURCES = $(zrank_mpi_SOURCES)
zrank_mpi_omp_CPPFLAGS= $(zrank_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
zrank_mpi_omp_CXXFLAGS= $(zrank_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
zrank_mpi_omp_LDFLAGS= $(zrank_mpi_LDFLAGS) # openmp?
zrank_mpi_omp_LDADD = $(zrank_mpi_LDADD) # openmp?


# unit tests
check_PROGRAMS = test-sparserow test-denserow test-rank-int test-rank-double

test_sparserow_SOURCES = test/test_sparserow.cpp
test_sparserow_CXXFLAGS = $(BOOST_CPPFLAGS)
test_sparserow_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_sparserow_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

test_denserow_SOURCES = test/test_denserow.cpp
test_denserow_CXXFLAGS = $(BOOST_CPPFLAGS)
test_denserow_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_denserow_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

test_rank_int_SOURCES = test/test_rank_int.cpp
test_rank_int_CXXFLAGS = $(BOOST_CPPFLAGS)
test_rank_int_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rank_int_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

test_rank_double_SOURCES = test/test_rank_double.cpp
test_rank_double_CXXFLAGS = $(BOOST_CPPFLAGS)
test_rank_double_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rank_double_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS = $(check_PROGRAMS)
