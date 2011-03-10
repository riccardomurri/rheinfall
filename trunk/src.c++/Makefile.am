## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rheinfall.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/types.hpp

# build example programs
bin_PROGRAMS = rank-int rank-double rank-mod 
EXTRA_PROGRAMS = rank-xint
if HAVE_BOOST_XINT
  bin_PROGRAMS += rank-xint
endif

EXTRA_PROGRAMS += rank-int32 rank-int64
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32
endif
if HAVE_INT64_T
  bin_PROGRAMS += rank-int64 rank-mod-int64
endif

EXTRA_PROGRAMS += rank-mpq rank-mpz
if HAVE_GMPXX
  bin_PROGRAMS += rank-mpq rank-mpz
endif

EXTRA_PROGRAMS += rank-int-omp rank-double-omp rank-mod-omp rank-mpq-omp rank-mpz-omp
EXTRA_PROGRAMS += rank-int32-omp 
EXTRA_PROGRAMS += rank-int64-omp rank-mod-int64-omp
if HAVE_OPENMP_CXX
  bin_PROGRAMS +=  rank-int-omp rank-mod-omp 
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32-omp
endif
if HAVE_INT64_T
  bin_PROGRAMS += rank-int64-omp rank-mod-int64-omp
endif
if HAVE_GMPXX 
  bin_PROGRAMS += rank-mpq-omp rank-mpz-omp
endif
endif
EXTRA_PROGRAMS += rank-int-mpi rank-double rank-mod-mpi rank-mpq-mpi rank-mpz-mpi
EXTRA_PROGRAMS += rank-int-mpi-omp rank-mod-mpi-omp  rank-mpq-mpi-omp rank-mpz-mpi-omp
EXTRA_PROGRAMS += rank-int32-mpi rank-int32-mpi-omp
EXTRA_PROGRAMS += rank-int64-mpi rank-int64-mpi-omp
if HAVE_MPICXX
  bin_PROGRAMS += rank-int-mpi rank-double-mpi rank-mod-mpi #rank-mpq-mpi rank-mpz-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += rank-int-mpi-omp rank-double-mpi-omp rank-mod-mpi-omp #rank-mpq-mpi-omp rank-mpz-mpi-omp
endif
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32-mpi
if HAVE_OPENMP_CXX
    bin_PROGRAMS += rank-int32-mpi-omp
endif
endif
endif


AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_srcdir)/src.c++/rheinfall $(BOOST_CPPFLAGS)
if WITH_GE
  AM_CPPFLAGS += -DWITH_GE
endif
AM_CXXFLAGS = 

# irank - example program using `intmax_t`/`long long` (whichever is wider) integer values
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
rank_int_mpi_omp_CPPFLAGS= $(rank_int_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int_mpi_omp_CXXFLAGS= $(rank_int_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int_mpi_omp_LDFLAGS= $(rank_int_mpi_LDFLAGS) # openmp?
rank_int_mpi_omp_LDADD = $(rank_int_mpi_LDADD) # openmp?

rank_int32_SOURCES = examples/rank.cpp
rank_int32_CPPFLAGS= $(AM_CPPFLAGS) -DWITH_INT32_VALUES
rank_int32_CXXFLAGS= $(AM_CXXFLAGS)
rank_int32_LDFLAGS= $(AM_LDFLAGS) 
rank_int32_LDADD = $(AM_LIBS) 

rank_int32_omp_SOURCES = $(rank_int32_SOURCES)
rank_int32_omp_CPPFLAGS = $(rank_int32_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int32_omp_CXXFLAGS = $(rank_int32_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int32_omp_LDFLAGS = $(rank_int32_LDFLAGS) # openmp?
rank_int32_omp_LDADD = $(rank_int32_LDADD) # openmp?

rank_int32_mpi_SOURCES = $(rank_int32_SOURCES)
rank_int32_mpi_CPPFLAGS = $(rank_int32_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int32_mpi_CXXFLAGS = $(rank_int32_CXXFLAGS) 
rank_int32_mpi_LDFLAGS = $(rank_int32_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int32_mpi_LDADD = $(rank_int32_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int32_mpi_omp_SOURCES = $(rank_int32_mpi_SOURCES)
rank_int32_mpi_omp_CPPFLAGS= $(rank_int32_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int32_mpi_omp_CXXFLAGS= $(rank_int32_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int32_mpi_omp_LDFLAGS= $(rank_int32_mpi_LDFLAGS) # openmp?
rank_int32_mpi_omp_LDADD = $(rank_int32_mpi_LDADD) # openmp?

rank_int64_SOURCES = examples/rank.cpp
rank_int64_CPPFLAGS= $(AM_CPPFLAGS) -DWITH_INT64_VALUES
rank_int64_CXXFLAGS= $(AM_CXXFLAGS)
rank_int64_LDFLAGS= $(AM_LDFLAGS) 
rank_int64_LDADD = $(AM_LIBS) 

rank_int64_omp_SOURCES = $(rank_int64_SOURCES)
rank_int64_omp_CPPFLAGS = $(rank_int64_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int64_omp_CXXFLAGS = $(rank_int64_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int64_omp_LDFLAGS = $(rank_int64_LDFLAGS) # openmp?
rank_int64_omp_LDADD = $(rank_int64_LDADD) # openmp?

rank_int64_mpi_SOURCES = $(rank_int64_SOURCES)
rank_int64_mpi_CPPFLAGS = $(rank_int64_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int64_mpi_CXXFLAGS = $(rank_int64_CXXFLAGS) 
rank_int64_mpi_LDFLAGS = $(rank_int64_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int64_mpi_LDADD = $(rank_int64_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int64_mpi_omp_SOURCES = $(rank_int64_mpi_SOURCES)
rank_int64_mpi_omp_CPPFLAGS= $(rank_int64_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_int64_mpi_omp_CXXFLAGS= $(rank_int64_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_int64_mpi_omp_LDFLAGS= $(rank_int64_mpi_LDFLAGS) # openmp?
rank_int64_mpi_omp_LDADD = $(rank_int64_mpi_LDADD) # openmp?

# rank-mod - example program using mod p values (over `long long` or `long`, whichever is wider)
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
rank_mod_mpi_omp_CPPFLAGS= $(rank_mod_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_mod_mpi_omp_CXXFLAGS= $(rank_mod_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mod_mpi_omp_LDFLAGS= $(rank_mod_mpi_LDFLAGS) # openmp?
rank_mod_mpi_omp_LDADD = $(rank_mod_mpi_LDADD) # openmp?

# drank - example program using `double`/`long double` values (whichever is wider)
rank_double_SOURCES = examples/rank.cpp
rank_double_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_DOUBLE_VALUES
rank_double_CXXFLAGS= $(AM_CXXFLAGS)
rank_double_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
rank_double_LDADD= $(AM_LIBS) $(GMPXX_LIBS)

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
rank_double_mpi_omp_CPPFLAGS= $(rank_double_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
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
rank_mpq_mpi_omp_CPPFLAGS= $(rank_mpq_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
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
rank_mpz_mpi_omp_CPPFLAGS= $(rank_mpz_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(OPENMP_CPPFLAGS)
rank_mpz_mpi_omp_CXXFLAGS= $(rank_mpz_mpi_CXXFLAGS) $(OPENMP_CXXFLAGS)
rank_mpz_mpi_omp_LDFLAGS= $(rank_mpz_mpi_LDFLAGS) # openmp?
rank_mpz_mpi_omp_LDADD = $(rank_mpz_mpi_LDADD) # openmp?

# example program using Boost.XInt infinite-precision integer values
rank_xint_SOURCES = examples/rank.cpp
rank_xint_CPPFLAGS= $(AM_CPPFLAGS)  -DWITH_XINT_VALUES
rank_xint_CXXFLAGS= $(AM_CXXFLAGS)
rank_xint_LDFLAGS= $(AM_LDFLAGS) $(BOOST_XINT_LDFLAGS)
rank_xint_LDADD = $(AM_LIBS) $(BOOST_XINT_LIBS)



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

# test matrices are downloaded from the web
data_DATA = M0,5-D4.sms M0,5-D5.sms M0,5-D6.sms M0,5-D7.sms M0,5-D8.sms \
		M0,6-D5.sms M0,6-D6.sms M0,6-D7.sms M0,6-D8.sms M0,6-D11.sms
$(data_DATA):
	@WGET@ -nv -nc -N \
		http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/Matrices/Mgn/$@.gz
	@GZIP@ -d $(srcdir)/$@.gz

distclean-local:
	rm -f $(data_DATA)
