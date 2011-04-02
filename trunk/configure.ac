## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

AC_PREREQ([2.65])
AC_INIT([rheinfall],[0.15.10],[riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c++/rheinfall/rheinfall.hpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AM_INIT_AUTOMAKE([1.10 no-define])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
gl_EARLY
AC_PROG_CXX
AC_PROG_INSTALL

AC_PATH_PROG([GZIP],  [gzip], [:])
AC_PATH_PROG([WGET],  [wget], [:])


# initialize GNU LibTool
LT_INIT


# Initialize GNUlib
AC_LANG([C])
gl_INIT


# check for OpenMP
AC_LANG_PUSH([C])
AC_OPENMP
AC_LANG_POP([C])
AM_CONDITIONAL([HAVE_OPENMP_CC], [test "_$OPENMP_CFLAGS" != '_'])

AC_LANG_PUSH([C++])
AC_OPENMP
AC_LANG_POP([C++])
AM_CONDITIONAL([HAVE_OPENMP_CXX], [test "_$OPENMP_CXXFLAGS" != '_'])


# check for MPI
AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi=no|yes|serialized|multiple],
            [Build MPI-enabled example programs. 
                If argument is 'yes' (default), the MPI C/C++ compilers must be
                reachable with the standard executable search path PATH.
                Argument values 'serialized' and 'multiple' behave as
                'yes' and additionally force use of MPI with threads,
                setting the threading level to MPI_THREAD_SERIALIZED
                or MPI_THREAD_MULTIPLE respectively.
                If the argument is 'no', then example programs that make 
                use of MPI will not be built.])],
    [],
    [with_mpi=yes])
if test "_$with_mpi" != '_no'; then :
    AC_LANG_PUSH([C])
    LX_FIND_MPI()
    AC_LANG_POP([C])
    AC_LANG_PUSH([C++])
    LX_FIND_MPI()
    AC_LANG_POP([C++])
fi
case "_$with_mpi" in #(
     _s*) AC_SUBST([RF_MPI_CPPFLAGS], [-DWITH_MPI_SERIALIZED]) ;; #(
     _m*) AC_SUBST([RF_MPI_CPPFLAGS], [-DWITH_MPI_MULTIPLE]) ;; #(
     *)   AC_SUBST([RF_MPI_CPPFLAGS], []) ;;
esac


# check for needed Boost libraries 
BOOST_REQUIRE([1.45]) # issend patch only tested starting 1.45.0
BOOST_STATIC
BOOST_SERIALIZATION
BOOST_STATICASSERT
BOOST_TEST([], 
    [have_boost_test=no;
     AC_MSG_WARN([Could not find a usable Boost.Test library; not compiling test programs.])],
    [have_boost_test=yes])
AM_CONDITIONAL([HAVE_BOOST_TEST], [test "_$have_boost_test" = '_yes'])
BOOST_XINT([AC_MSG_WARN([cannot find Boost.XInt; no XInt test programs will be compiled.])])
AM_CONDITIONAL([HAVE_BOOST_XINT], 
               [test "_$ac_cv_header_boost_xint_integer_hpp" = '_yes'])

if test "_$have_CXX_mpi" = '_yes'; then
    BOOST_MPI([], 
        [# disable MPI/C++ altogether
           have_CXX_mpi=no;
           AC_MSG_WARN([Could not find a usable Boost.MPI library; not compiling MPI example C++ programs.])],
        [# check for `issend`
         AC_MSG_CHECKING([whether Boost.MPI provides the issend() function])
         AC_LANG_PUSH([C++])
         save_CXXFLAGS="$CXXFLAGS"
         CXXFLAGS="$CXXFLAGS $MPI_CXXFLAGS"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#include <boost/mpi.hpp>
int main(int argc, char** argv) 
{
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator world;
        int i = 0;
        world.issend(0, 0, i);
        return 0;
}
                                            ])],
                           [AC_MSG_RESULT([yes])],
                           [# disable MPI/C++ altogether
                            AC_MSG_RESULT([no])
                            have_CXX_mpi=no;
                            AC_MSG_WARN([Boost.MPI does not provide the issend() function; not compiling MPI example C++ programs.])])
         CXXFLAGS="$save_CXXFLAGS"
         AC_LANG_POP([C++])
        ], 
        [$MPI_CXXFLAGS $BOOST_CPPFLAGS],
        [$MPI_CXXFLAGS $BOOST_CPPFLAGS],
        [$MPI_CXXLDFLAGS $BOOST_SERIALIZATION_LDFLAGS],
        [$MPI_CXXLDFLAGS $BOOST_SERIALIZATION_LIBS])
fi

AM_CONDITIONAL([HAVE_MPICC], [test "_$have_C_mpi" = '_yes'])
AM_CONDITIONAL([HAVE_MPICXX], [test "_$have_CXX_mpi" = '_yes'])


# Check for GMP/GMPXX
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=yes|no|@<:@path@:>@],
            [Build example programs requiring the GMP library. 
                If argument is 'yes' or empty, the library must be 
                reachable with the standard compiler search path; 
                otherwise the GMP include files have to be in 
                <path>/include and the compiled library has to 
                reside in <path>/lib.  If the argument is 'no', then
                example programs that make use of GMP will not be built.])],
    [ # option '--with-gmp' given
      if test "_$withval" = _yes ; then :
          # require GMP and search for it in default system paths
          AX_CHECK_GMP
          AX_CHECK_GMPXX
      elif test "_$withval" != _no ; then :
          AX_CHECK_GMP([], "$withval")
          AX_CHECK_GMPXX
      fi
    ])
AM_CONDITIONAL([HAVE_GMPXX], [test "_$have_gmpxx" = '_yes'])


# Check for Grid Engine (Sun/Oracle)
AC_ARG_WITH([ge],
    [AS_HELP_STRING([--with-ge=yes|no],
            [Print abort message and exit upon receiving "kill" notification from GE.])],
    [],
    [with_ge=yes])
AM_CONDITIONAL([WITH_GE], [test "_$with_ge" = '_yes'])


# Checks for header files.
AC_CHECK_HEADERS([inttypes.h malloc.h stdint.h stdlib.h string.h sys/time.h sys/types.h unistd.h])

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([ext/malloc_allocator.h])
AC_LANG_POP([C++])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AX_CXX_RESTRICT_THIS
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_INTMAX_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_LONG_DOUBLE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AM_CONDITIONAL([HAVE_INT32_T], [test "_$ac_cv_c_int32_t" = '_yes'])
AM_CONDITIONAL([HAVE_INT64_T], [test "_$ac_cv_c_int64_t" = '_yes'])
AM_CONDITIONAL([HAVE_LONG_LONG_INT], [test "_$ac_cv_type_long_long_int" = '_yes'])
AM_CONDITIONAL([HAVE_LONG_DOUBLE], [test "_$ac_cv_type_long_double" = '_yes'])


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([backtrace backtrace_symbols_fd gettimeofday memmove])


# Final output
AC_CONFIG_FILES([
        src.c++/Makefile
        src.c/Makefile
        src.c/gnulib/Makefile
])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
