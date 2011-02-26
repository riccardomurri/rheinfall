## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

AC_PREREQ([2.63])
AC_INIT([rheinfall], [0.10.15], [riccardo.murri@gmail.com])
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


# check for OpenMP
AC_LANG_PUSH([C])
AC_OPENMP
AC_LANG_POP([C])
AM_CONDITIONAL([HAVE_OPENMP_CC], [test "x$OPENMP_CFLAGS" != x])

AC_LANG_PUSH([C++])
AC_OPENMP
AC_LANG_POP([C++])
AM_CONDITIONAL([HAVE_OPENMP_CXX], [test "x$OPENMP_CXXFLAGS" != x])


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
AM_CONDITIONAL([HAVE_MPICC], [test "x$have_C_mpi" != no])
AM_CONDITIONAL([HAVE_MPICXX], [test "x$have_CXX_mpi" != no])


# check for needed Boost libraries
BOOST_REQUIRE([1.35]) # Boost.MPI released with Boost 1.35.0
BOOST_STATIC
BOOST_SERIALIZATION
BOOST_MPI([], 
    [], 
    [# disable MPI/C++ altogether
        have_CXX_mpi=no;
        AC_MSG_WARN([Could not find a usable Boost.MPI library; not compiling MPI example C++ programs.])],
    [$MPI_CXXFLAGS],
    [$MPI_CXXFLAGS],
    [$MPI_CXXLDFLAGS],
    [$MPI_CXXLDFLAGS])
BOOST_STATICASSERT
BOOST_TEST([], 
           [have_boost_test=yes], 
           [have_boost_test=no;
            AC_MSG_WARN([Could not find a usable Boost.Test library; not compiling test programs.])])
AM_CONDITIONAL([HAVE_BOOST_TEST], [test "x$have_boost_test" != 'xno'])


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
      if test "$withval" = yes ; then :
          # require GMP and search for it in default system paths
          AX_CHECK_GMP
          AX_CHECK_GMPXX
      elif test "$withval" != no ; then :
          AX_CHECK_GMP([], "$withval")
          AX_CHECK_GMPXX
      fi
    ])
AM_CONDITIONAL([HAVE_GMPXX], [test "x$have_gmpxx" != 'xno'])


# check for OpenMP
AC_OPENMP

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_INTMAX_T


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
