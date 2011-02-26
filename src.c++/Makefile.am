## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rheinfall.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/mathdefs.hpp

# build example programs
bin_PROGRAMS = rank-int rank-double rank-mpq rank-mpz
EXTRA_PROGRAMS = rank-int-omp rank-double-omp  rank-mpq-omp rank-mpz-omp
if HAVE_OPENMP_CXX
  bin_PROGRAMS +=  rank-int-omp rank-double-omp  rank-mpq-omp rank-mpz-omp
endif
EXTRA_PROGRAMS += rank-int-mpi rank-double-mpi  rank-mpq-mpi rank-mpz-mpi
EXTRA_PROGRAMS += rank-int-mpi-omp rank-double-mpi-omp  rank-mpq-mpi-omp rank-mpz-mpi-omp
if HAVE_MPICXX
  bin_PROGRAMS += rank-int-mpi rank-double-mpi  #rank-mpq-mpi rank-mpz-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += rank-int-mpi-omp rank-double-mpi-omp  #rank-mpq-mpi-omp rank-mpz-mpi-omp
endif
endif

if HAVE_OPENMP_CXX
omp: rank-int-omp rank-double-omp  rank-mpq-omp rank-mpz-omp
endif
if HAVE_MPICXX
mpi: rank-int-mpi rank-double-mpi  #rank-mpq-mpi rank-mpz-mpi
if HAVE_OPENMP_CXX
mpi-omp: rank-int-mpi-omp rank-double-mpi-omp  #rank-mpq-mpi-omp rank-mpz-mpi-omp
endif
endif


AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c++/rheinfall $(BOOST_CPPFLAGS)
AM_CXXFLAGS = 

# irank - example program using `intmax_t`/`long long` integer values
rank_int_SOURCES = examples/rank.cpp
rank_int_CPPFLAGS= $(AM_CPPFLAGS) -DWITH_INT_VALUES
rank_int_CXXFLAGS= $(AM_CXXFLAGS)
rank_int_LDFLAGS= $(AM_LDFLAGS) 
rank_int_LDADD = $(AM_LIBS) 

rank_int_omp_SOURCES = $(rank_int_SOURCES)
rank_int_omp_CPPFLAGS = $(rank_int_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int_omp_CXXFLAGS = $(rank_int_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int_omp_LDFLAGS = $(rank_int_LDFLAGS) # openmp?
rank_int_omp_LDADD = $(rank_int_LDADD) # openmp?

rank_int_mpi_SOURCES = $(rank_int_SOURCES)
rank_int_mpi_CPPFLAGS = $(rank_int_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int_mpi_CXXFLAGS = $(rank_int_CXXFLAGS) 
rank_int_mpi_LDFLAGS = $(rank_int_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int_mpi_LDADD = $(rank_int_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int_mpi_omp_SOURCES = $(rank_int_mpi_SOURCES)
rank_int_mpi_omp_CPPFLAGS= $(rank_int_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
rank_int_mpi_omp_CXXFLAGS= $(rank_int_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int_mpi_omp_LDFLAGS= $(rank_int_mpi_LDFLAGS) # openmp?
rank_int_mpi_omp_LDADD = $(rank_int_mpi_LDADD) # openmp?

# drank - example program using `long double`/`double` floating-point values
rank_double_SOURCES = examples/rank.cpp
rank_double_CPPFLAGS= $(AM_CPPFLAGS)  -DWITH_DOUBLE_VALUES
rank_double_CXXFLAGS= $(AM_CXXFLAGS)
rank_double_LDFLAGS= $(AM_LDFLAGS) 
rank_double_LDADD = $(AM_LIBS) 

rank_double_omp_SOURCES = $(rank_double_SOURCES)
rank_double_omp_CPPFLAGS = $(rank_double_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_double_omp_CXXFLAGS = $(rank_double_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_double_omp_LDFLAGS = $(rank_double_LDFLAGS) # openmp?
rank_double_omp_LDADD = $(rank_double_LDADD) # openmp?

rank_double_mpi_SOURCES = $(rank_double_SOURCES)
rank_double_mpi_CPPFLAGS = $(rank_double_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_double_mpi_CXXFLAGS = $(rank_double_CXXFLAGS)
rank_double_mpi_LDFLAGS = $(rank_double_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_double_mpi_LDADD = $(rank_double_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_double_mpi_omp_SOURCES = $(rank_double_mpi_SOURCES)
rank_double_mpi_omp_CPPFLAGS= $(rank_double_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
rank_double_mpi_omp_CXXFLAGS= $(rank_double_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_double_mpi_omp_LDFLAGS= $(rank_double_mpi_LDFLAGS) # openmp?
rank_double_mpi_omp_LDADD = $(rank_double_mpi_LDADD) # openmp?

# qrank - example program using GMP `mpq` infinite precision rational values
rank_mpq_SOURCES = examples/rank.cpp
rank_mpq_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_MPQ_VALUES
rank_mpq_CXXFLAGS= $(AM_CXXFLAGS)
rank_mpq_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
rank_mpq_LDADD= $(AM_LIBS) $(GMPXX_LIBS)

rank_mpq_omp_SOURCES = $(rank_mpq_SOURCES)
rank_mpq_omp_CPPFLAGS = $(rank_mpq_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_mpq_omp_CXXFLAGS = $(rank_mpq_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mpq_omp_LDFLAGS = $(rank_mpq_LDFLAGS) # openmp?
rank_mpq_omp_LDADD = $(rank_mpq_LDADD) # openmp?

rank_mpq_mpi_SOURCES = $(rank_mpq_SOURCES)
rank_mpq_mpi_CPPFLAGS = $(rank_mpq_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mpq_mpi_CXXFLAGS = $(rank_mpq_CXXFLAGS)
rank_mpq_mpi_LDFLAGS = $(rank_mpq_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mpq_mpi_LDADD = $(rank_mpq_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mpq_mpi_omp_SOURCES = $(rank_mpq_mpi_SOURCES)
rank_mpq_mpi_omp_CPPFLAGS= $(rank_mpq_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
rank_mpq_mpi_omp_CXXFLAGS= $(rank_mpq_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mpq_mpi_omp_LDFLAGS= $(rank_mpq_mpi_LDFLAGS) # openmp?
rank_mpq_mpi_omp_LDADD = $(rank_mpq_mpi_LDADD) # openmp?

# zrank - example program using GMP `mpz_t` infinite-precision integer values
rank_mpz_SOURCES = examples/rank.cpp
rank_mpz_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_MPZ_VALUES
rank_mpz_CXXFLAGS= $(AM_CXXFLAGS)
rank_mpz_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
rank_mpz_LDADD = $(AM_LIBS) $(GMPXX_LIBS)

rank_mpz_omp_SOURCES = $(rank_mpz_SOURCES)
rank_mpz_omp_CPPFLAGS = $(rank_mpz_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_mpz_omp_CXXFLAGS = $(rank_mpz_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mpz_omp_LDFLAGS = $(rank_mpz_LDFLAGS) # openmp?
rank_mpz_omp_LDADD = $(rank_mpz_LDADD) # openmp?

rank_mpz_mpi_SOURCES = $(rank_mpz_SOURCES)
rank_mpz_mpi_CPPFLAGS = $(rank_mpz_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mpz_mpi_CXXFLAGS = $(rank_mpz_CXXFLAGS)
rank_mpz_mpi_LDFLAGS = $(rank_mpz_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mpz_mpi_LDADD = $(rank_mpz_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mpz_mpi_omp_SOURCES = $(rank_mpz_mpi_SOURCES)
rank_mpz_mpi_omp_CPPFLAGS= $(rank_mpz_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
rank_mpz_mpi_omp_CXXFLAGS= $(rank_mpz_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mpz_mpi_omp_LDFLAGS= $(rank_mpz_mpi_LDFLAGS) # openmp?
rank_mpz_mpi_omp_LDADD = $(rank_mpz_mpi_LDADD) # openmp?


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
