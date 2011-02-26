## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rheinfall.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/types.hpp

# build example programs
bin_PROGRAMS = rank-int rank-mod rank-mpq rank-mpz
EXTRA_PROGRAMS = rank-int-omp rank-mod-omp  rank-mpq-omp rank-mpz-omp
if HAVE_OPENMP_CXX
  bin_PROGRAMS +=  rank-int-omp rank-mod-omp  rank-mpq-omp rank-mpz-omp
endif
EXTRA_PROGRAMS += rank-int-mpi rank-mod-mpi  rank-mpq-mpi rank-mpz-mpi
EXTRA_PROGRAMS += rank-int-mpi-omp rank-mod-mpi-omp  rank-mpq-mpi-omp rank-mpz-mpi-omp
if HAVE_MPICXX
  bin_PROGRAMS += rank-int-mpi #rank-mod-mpi rank-mpq-mpi rank-mpz-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += rank-int-mpi-omp #rank-mod-mpi-omp rank-mpq-mpi-omp rank-mpz-mpi-omp
endif
endif

if HAVE_OPENMP_CXX
omp: rank-int-omp rank-mod-omp  rank-mpq-omp rank-mpz-omp
endif
if HAVE_MPICXX
mpi: rank-int-mpi #rank-mod-mpi rank-mpq-mpi rank-mpz-mpi
if HAVE_OPENMP_CXX
mpi-omp: rank-int-mpi-omp #rank-mod-mpi-omp rank-mpq-mpi-omp rank-mpz-mpi-omp
endif
endif


AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c++/rheinfall $(BOOST_CPPFLAGS)
if WITH_GE
  AM_CPPFLAGS += -DWITH_GE
endif
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

# rank-mod - example program using mod p values (over `long long` or `long`)
rank_mod_SOURCES = examples/rank.cpp
rank_mod_CPPFLAGS= $(AM_CPPFLAGS)  -DWITH_MODULAR_VALUES
rank_mod_CXXFLAGS= $(AM_CXXFLAGS)
rank_mod_LDFLAGS= $(AM_LDFLAGS) 
rank_mod_LDADD = $(AM_LIBS) 

rank_mod_omp_SOURCES = $(rank_mod_SOURCES)
rank_mod_omp_CPPFLAGS = $(rank_mod_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_mod_omp_CXXFLAGS = $(rank_mod_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mod_omp_LDFLAGS = $(rank_mod_LDFLAGS) # openmp?
rank_mod_omp_LDADD = $(rank_mod_LDADD) # openmp?

rank_mod_mpi_SOURCES = $(rank_mod_SOURCES)
rank_mod_mpi_CPPFLAGS = $(rank_mod_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mod_mpi_CXXFLAGS = $(rank_mod_CXXFLAGS)
rank_mod_mpi_LDFLAGS = $(rank_mod_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mod_mpi_LDADD = $(rank_mod_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mod_mpi_omp_SOURCES = $(rank_mod_mpi_SOURCES)
rank_mod_mpi_omp_CPPFLAGS= $(rank_mod_mpi_CPPFLAGS) -DWITH_MPI_SERIALIZED $(OPENMP_CPPFLAGS)
rank_mod_mpi_omp_CXXFLAGS= $(rank_mod_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mod_mpi_omp_LDFLAGS= $(rank_mod_mpi_LDFLAGS) # openmp?
rank_mod_mpi_omp_LDADD = $(rank_mod_mpi_LDADD) # openmp?

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
check_PROGRAMS = test-sparserow test-denserow test-rank-int test-rank-mod

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

test_rank_mod_SOURCES = test/test_rank_mod.cpp
test_rank_mod_CXXFLAGS = $(BOOST_CPPFLAGS)
test_rank_mod_LDFLAGS= $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
test_rank_mod_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIBS)

TESTS = $(check_PROGRAMS)
