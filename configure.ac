## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

AC_PREREQ([2.63])
AC_INIT([rheinfall], [0.11.26], [riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c++/rheinfall/rheinfall.hpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AM_INIT_AUTOMAKE([1.10 no-define])
AM_CONFIG_HEADER(config.h)


# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
gl_EARLY
AC_PROG_CXX
AC_PROG_INSTALL


# initialize GNU LibTool
LT_INIT


# Initialize GNUlib
AC_LANG([C])
gl_INIT


# DON'T check for OpenMP -- the code won't run... (to be fixed someday).
#AC_LANG_PUSH([C])
#AC_OPENMP
#AC_LANG_POP([C])
AM_CONDITIONAL([HAVE_OPENMP_CC], [test "_$OPENMP_CFLAGS" != '_'])

#AC_LANG_PUSH([C++])
#AC_OPENMP
#AC_LANG_POP([C++])
AM_CONDITIONAL([HAVE_OPENMP_CXX], [test "_$OPENMP_CXXFLAGS" != '_'])


# check for MPI
AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi=yes|no],
            [Build MPI-enabled example programs. 
                If argument is 'yes' (default), the MPI C/C++ compilers must be
                reachable with the standard executable search path PATH.
                If the argument is 'no', then example programs that make 
                use of MPI will not be built.])],
    [],
    [with_mpi=yes])
if test "x$with_mpi" = 'xyes'; then :
    AC_LANG_PUSH([C])
    LX_FIND_MPI()
    AC_LANG_POP([C])
    AC_LANG_PUSH([C++])
    LX_FIND_MPI()
    AC_LANG_POP([C++])
fi
AM_CONDITIONAL([HAVE_MPICC], [test "_$have_C_mpi" = '_yes'])
AM_CONDITIONAL([HAVE_MPICXX], [test "_$have_CXX_mpi" = '_yes'])


# check for needed Boost libraries 
BOOST_REQUIRE([1.35]) # Boost.MPI released with Boost 1.35.0
BOOST_STATIC
BOOST_SERIALIZATION
BOOST_TEST([], 
    [have_boost_test=no;
        AC_MSG_WARN([Could not find a usable Boost.Test library; not compiling test programs.])],
    [have_boost_test=yes])
AM_CONDITIONAL([HAVE_BOOST_TEST], [test "_$have_boost_test" = '_yes'])
if test "_$have_CXX_mpi" = '_yes'; then
    if test "_$boost_major_version" = "_144"; then
        # see: http://lists.boost.org/boost-users/2010/11/63980.php
        AC_MSG_ERROR([There is a known issue with Boost.MPI 1.44, please use 1.43 or a later version.])
    fi
    BOOST_MPI([], 
        [# disable MPI/C++ altogether
            have_CXX_mpi=no;
            AC_MSG_WARN([Could not find a usable Boost.MPI library; not compiling MPI example C++ programs.])],
        [], 
        [$MPI_CXXFLAGS],
        [$MPI_CXXFLAGS],
        [$MPI_CXXLDFLAGS],
        [$MPI_CXXLDFLAGS])
fi
BOOST_STATICASSERT


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
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AX_CXX_RESTRICT_THIS
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_INTMAX_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_LONG_DOUBLE
AC_CHECK_TYPES([int128_t], [],
    [AC_CHECK_TYPES([__int128_t],
            [AC_DEFINE([HAVE_INT128_T])
             AC_DEFINE([int128_t], [__int128_t], [Name of the 128-bit integer type.])
             ])])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove])


# Final output
AC_CONFIG_FILES([
        src.c++/Makefile
        src.c/Makefile
        src.c/gnulib/Makefile
])
AC_OUTPUT(Makefile)
