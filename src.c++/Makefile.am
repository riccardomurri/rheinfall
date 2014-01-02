## Process this file with automake to produce Makefile.in

# this is an include-only library
nobase_include_HEADERS = rheinfall/rank.hpp rheinfall/row.hpp rheinfall/sparserow.hpp rheinfall/denserow.hpp rheinfall/types.hpp rheinfall/types/gmpxx.hpp rheinfall/types/modular.hpp

# build example programs
bin_PROGRAMS = rank-int rank-mod lu

EXTRA_PROGRAMS = rank-gfq
if HAVE_GIVARO
  bin_PROGRAMS += rank-gfq
endif

EXTRA_PROGRAMS += rank-int32 rank-int64
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32
endif
if HAVE_INT64_T
  bin_PROGRAMS += rank-int64
endif

EXTRA_PROGRAMS += rank-mpq rank-mpz
if HAVE_GMPXX
  bin_PROGRAMS += rank-mpq rank-mpz
endif

EXTRA_PROGRAMS += rank-int-tbb rank-mod-tbb rank-mpq-tbb rank-mpz-tbb
EXTRA_PROGRAMS += rank-int32-tbb 
EXTRA_PROGRAMS += rank-int64-tbb
if HAVE_TBB
  bin_PROGRAMS +=  rank-int-tbb rank-mod-tbb lu-tbb
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32-tbb
endif
if HAVE_INT64_T
  bin_PROGRAMS += rank-int64-tbb
endif
if HAVE_GMPXX 
  bin_PROGRAMS += rank-mpq-tbb rank-mpz-tbb
endif
endif
EXTRA_PROGRAMS += rank-int-mpi rank-mod-mpi rank-mpq-mpi rank-mpz-mpi
EXTRA_PROGRAMS += rank-int-mpi-tbb rank-mod-mpi-tbb  rank-mpq-mpi-tbb rank-mpz-mpi-tbb
EXTRA_PROGRAMS += rank-int32-mpi rank-int32-mpi-tbb
EXTRA_PROGRAMS += rank-int64-mpi rank-int64-mpi-tbb
if HAVE_MPICXX
  bin_PROGRAMS += rank-int-mpi rank-mod-mpi lu-mpi #rank-mpq-mpi rank-mpz-mpi
if HAVE_TBB
    bin_PROGRAMS += rank-int-mpi-tbb rank-mod-mpi-tbb lu-mpi-tbb #rank-mpq-mpi-tbb rank-mpz-mpi-tbb
endif
if HAVE_INT32_T
  bin_PROGRAMS += rank-int32-mpi
if HAVE_TBB
    bin_PROGRAMS += rank-int32-mpi-tbb
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

rank_gfq_SOURCES = examples/rank.cpp
rank_gfq_CPPFLAGS= $(AM_CPPFLAGS) $(GIVARO_CPPFLAGS) -DWITH_GFQ_VALUES -I/home/rmurri/my/linbox/linbox-1.3.2/sw/include
rank_gfq_CXXFLAGS= $(AM_CXXFLAGS) $(GIVARO_CXXFLAGS)
rank_gfq_LDFLAGS= $(AM_LDFLAGS) $(GIVARO_LDFLAGS)
rank_gfq_LDADD = $(AM_LIBS) $(GIVARO_LIBS)

rank_int_tbb_SOURCES = $(rank_int_SOURCES)
rank_int_tbb_CPPFLAGS = $(rank_int_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int_tbb_CXXFLAGS = $(rank_int_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int_tbb_LDFLAGS = $(rank_int_LDFLAGS) $(TBB_LDFLAGS)
rank_int_tbb_LDADD = $(rank_int_LDADD) $(TBB_LIBS)

rank_int_mpi_SOURCES = $(rank_int_SOURCES)
rank_int_mpi_CPPFLAGS = $(rank_int_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int_mpi_CXXFLAGS = $(rank_int_CXXFLAGS) 
rank_int_mpi_LDFLAGS = $(rank_int_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int_mpi_LDADD = $(rank_int_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int_mpi_tbb_SOURCES = $(rank_int_mpi_SOURCES)
rank_int_mpi_tbb_CPPFLAGS= $(rank_int_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int_mpi_tbb_CXXFLAGS= $(rank_int_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int_mpi_tbb_LDFLAGS= $(rank_int_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_int_mpi_tbb_LDADD = $(rank_int_mpi_LDADD) $(TBB_LIBS)

rank_int32_SOURCES = examples/rank.cpp
rank_int32_CPPFLAGS= $(AM_CPPFLAGS) -DWITH_INT32_VALUES
rank_int32_CXXFLAGS= $(AM_CXXFLAGS)
rank_int32_LDFLAGS= $(AM_LDFLAGS) 
rank_int32_LDADD = $(AM_LIBS) 

rank_int32_tbb_SOURCES = $(rank_int32_SOURCES)
rank_int32_tbb_CPPFLAGS = $(rank_int32_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int32_tbb_CXXFLAGS = $(rank_int32_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int32_tbb_LDFLAGS = $(rank_int32_LDFLAGS) $(TBB_LDFLAGS)
rank_int32_tbb_LDADD = $(rank_int32_LDADD) $(TBB_LIBS)

rank_int32_mpi_SOURCES = $(rank_int32_SOURCES)
rank_int32_mpi_CPPFLAGS = $(rank_int32_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int32_mpi_CXXFLAGS = $(rank_int32_CXXFLAGS) 
rank_int32_mpi_LDFLAGS = $(rank_int32_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int32_mpi_LDADD = $(rank_int32_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int32_mpi_tbb_SOURCES = $(rank_int32_mpi_SOURCES)
rank_int32_mpi_tbb_CPPFLAGS= $(rank_int32_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int32_mpi_tbb_CXXFLAGS= $(rank_int32_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int32_mpi_tbb_LDFLAGS= $(rank_int32_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_int32_mpi_tbb_LDADD = $(rank_int32_mpi_LDADD) $(TBB_LIBS)

rank_int64_SOURCES = examples/rank.cpp
rank_int64_CPPFLAGS= $(AM_CPPFLAGS) -DWITH_INT64_VALUES
rank_int64_CXXFLAGS= $(AM_CXXFLAGS)
rank_int64_LDFLAGS= $(AM_LDFLAGS) 
rank_int64_LDADD = $(AM_LIBS) 

rank_int64_tbb_SOURCES = $(rank_int64_SOURCES)
rank_int64_tbb_CPPFLAGS = $(rank_int64_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int64_tbb_CXXFLAGS = $(rank_int64_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int64_tbb_LDFLAGS = $(rank_int64_LDFLAGS) $(TBB_LDFLAGS)
rank_int64_tbb_LDADD = $(rank_int64_LDADD) $(TBB_LIBS)

rank_int64_mpi_SOURCES = $(rank_int64_SOURCES)
rank_int64_mpi_CPPFLAGS = $(rank_int64_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_int64_mpi_CXXFLAGS = $(rank_int64_CXXFLAGS) 
rank_int64_mpi_LDFLAGS = $(rank_int64_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_int64_mpi_LDADD = $(rank_int64_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_int64_mpi_tbb_SOURCES = $(rank_int64_mpi_SOURCES)
rank_int64_mpi_tbb_CPPFLAGS= $(rank_int64_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_int64_mpi_tbb_CXXFLAGS= $(rank_int64_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_int64_mpi_tbb_LDFLAGS= $(rank_int64_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_int64_mpi_tbb_LDADD = $(rank_int64_mpi_LDADD) $(TBB_LIBS)

# rank-mod - example program using mod p values (over `long long` or `long`, whichever is wider)
rank_mod_SOURCES = examples/rank.cpp
rank_mod_CPPFLAGS= $(AM_CPPFLAGS)  -DWITH_MODULAR_VALUES
rank_mod_CXXFLAGS= $(AM_CXXFLAGS)
rank_mod_LDFLAGS= $(AM_LDFLAGS) 
rank_mod_LDADD = $(AM_LIBS) 

rank_mod_tbb_SOURCES = $(rank_mod_SOURCES)
rank_mod_tbb_CPPFLAGS = $(rank_mod_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mod_tbb_CXXFLAGS = $(rank_mod_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mod_tbb_LDFLAGS = $(rank_mod_LDFLAGS) $(TBB_LDFLAGS)
rank_mod_tbb_LDADD = $(rank_mod_LDADD) $(TBB_LIBS)

rank_mod_mpi_SOURCES = $(rank_mod_SOURCES)
rank_mod_mpi_CPPFLAGS = $(rank_mod_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mod_mpi_CXXFLAGS = $(rank_mod_CXXFLAGS)
rank_mod_mpi_LDFLAGS = $(rank_mod_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mod_mpi_LDADD = $(rank_mod_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mod_mpi_tbb_SOURCES = $(rank_mod_mpi_SOURCES)
rank_mod_mpi_tbb_CPPFLAGS= $(rank_mod_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mod_mpi_tbb_CXXFLAGS= $(rank_mod_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mod_mpi_tbb_LDFLAGS= $(rank_mod_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_mod_mpi_tbb_LDADD = $(rank_mod_mpi_LDADD) $(TBB_LIBS)

# qrank - example program using GMP `mpq` infinite precision rational values
rank_mpq_SOURCES = examples/rank.cpp
rank_mpq_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_MPQ_VALUES
rank_mpq_CXXFLAGS= $(AM_CXXFLAGS)
rank_mpq_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
rank_mpq_LDADD= $(AM_LIBS) $(GMPXX_LIBS)

rank_mpq_tbb_SOURCES = $(rank_mpq_SOURCES)
rank_mpq_tbb_CPPFLAGS = $(rank_mpq_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mpq_tbb_CXXFLAGS = $(rank_mpq_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mpq_tbb_LDFLAGS = $(rank_mpq_LDFLAGS) $(TBB_LDFLAGS)
rank_mpq_tbb_LDADD = $(rank_mpq_LDADD) $(TBB_LIBS)

rank_mpq_mpi_SOURCES = $(rank_mpq_SOURCES)
rank_mpq_mpi_CPPFLAGS = $(rank_mpq_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mpq_mpi_CXXFLAGS = $(rank_mpq_CXXFLAGS)
rank_mpq_mpi_LDFLAGS = $(rank_mpq_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mpq_mpi_LDADD = $(rank_mpq_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mpq_mpi_tbb_SOURCES = $(rank_mpq_mpi_SOURCES)
rank_mpq_mpi_tbb_CPPFLAGS= $(rank_mpq_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mpq_mpi_tbb_CXXFLAGS= $(rank_mpq_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mpq_mpi_tbb_LDFLAGS= $(rank_mpq_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_mpq_mpi_tbb_LDADD = $(rank_mpq_mpi_LDADD) $(TBB_LIBS)

# zrank - example program using GMP `mpz_t` infinite-precision integer values
rank_mpz_SOURCES = examples/rank.cpp
rank_mpz_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_MPZ_VALUES
rank_mpz_CXXFLAGS= $(AM_CXXFLAGS)
rank_mpz_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
rank_mpz_LDADD = $(AM_LIBS) $(GMPXX_LIBS)

rank_mpz_tbb_SOURCES = $(rank_mpz_SOURCES)
rank_mpz_tbb_CPPFLAGS = $(rank_mpz_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mpz_tbb_CXXFLAGS = $(rank_mpz_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mpz_tbb_LDFLAGS = $(rank_mpz_LDFLAGS) $(TBB_LDFLAGS)
rank_mpz_tbb_LDADD = $(rank_mpz_LDADD) $(TBB_LIBS)

rank_mpz_mpi_SOURCES = $(rank_mpz_SOURCES)
rank_mpz_mpi_CPPFLAGS = $(rank_mpz_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
rank_mpz_mpi_CXXFLAGS = $(rank_mpz_CXXFLAGS)
rank_mpz_mpi_LDFLAGS = $(rank_mpz_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
rank_mpz_mpi_LDADD = $(rank_mpz_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

rank_mpz_mpi_tbb_SOURCES = $(rank_mpz_mpi_SOURCES)
rank_mpz_mpi_tbb_CPPFLAGS= $(rank_mpz_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
rank_mpz_mpi_tbb_CXXFLAGS= $(rank_mpz_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
rank_mpz_mpi_tbb_LDFLAGS= $(rank_mpz_mpi_LDFLAGS) $(TBB_LDFLAGS)
rank_mpz_mpi_tbb_LDADD = $(rank_mpz_mpi_LDADD) $(TBB_LIBS)


# lu - example program using `double`/`long double` values (whichever is wider)
lu_SOURCES = examples/lu.cpp
lu_CPPFLAGS= $(AM_CPPFLAGS) $(GMPXX_CPPFLAGS) -DWITH_DOUBLE_VALUES
lu_CXXFLAGS= $(AM_CXXFLAGS)
lu_LDFLAGS= $(AM_LDFLAGS) $(GMPXX_LDFLAGS)
lu_LDADD= $(AM_LIBS) $(GMPXX_LIBS)

lu_tbb_SOURCES = $(lu_SOURCES)
lu_tbb_CPPFLAGS = $(lu_CPPFLAGS) $(TBB_CPPFLAGS)
lu_tbb_CXXFLAGS = $(lu_CXXFLAGS) $(TBB_CXXFLAGS)
lu_tbb_LDFLAGS = $(lu_LDFLAGS) $(TBB_LDFLAGS)
lu_tbb_LDADD = $(lu_LDADD) $(TBB_LIBS)

lu_mpi_SOURCES = $(lu_SOURCES)
lu_mpi_CPPFLAGS = $(lu_CPPFLAGS) -DWITH_MPI $(MPI_CXXFLAGS)
lu_mpi_CXXFLAGS = $(lu_CXXFLAGS)
lu_mpi_LDFLAGS = $(lu_LDFLAGS) $(BOOST_MPI_LDFLAGS) $(BOOST_SERIALIZATION_LDFLAGS) $(MPI_CXXLDFLAGS)
lu_mpi_LDADD = $(lu_LDADD) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LIBS) $(MPI_CXXLIBS)

lu_mpi_tbb_SOURCES = $(lu_mpi_SOURCES)
lu_mpi_tbb_CPPFLAGS= $(lu_mpi_CPPFLAGS) $(RF_MPI_CPPFLAGS) $(TBB_CPPFLAGS)
lu_mpi_tbb_CXXFLAGS= $(lu_mpi_CXXFLAGS) $(TBB_CXXFLAGS)
lu_mpi_tbb_LDFLAGS= $(lu_mpi_LDFLAGS) $(TBB_LDFLAGS)
lu_mpi_tbb_LDADD = $(lu_mpi_LDADD) $(TBB_LIBS)


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
	@WGET@ -nv -nc \
		http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/Matrices/Mgn/$@.gz
	@GZIP@ -d $(srcdir)/$@.gz

distclean-local:
	rm -f $(data_DATA)
