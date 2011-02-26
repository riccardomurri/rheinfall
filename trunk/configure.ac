## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

#AC_PREREQ([2.65])
AC_INIT([rheinfall], [0.10.15], [riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c++/rheinfall/rheinfall.hpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])


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
            [AS_HELP_STRING([--with-mpi],
              [prefix of MPI libraries @<:@guess@:>@])],
            [],
            [with_mpi=guess])
AM_CONDITIONAL([HAVE_MPICC], [test "x$have_mpicc" != no])
AM_CONDITIONAL([HAVE_MPICXX], [test "x$have_mpicxx" != no])
AS_IF([test "x$with_mpi" != xno],
  [# check for MPI compiler support
   save_CXX="$CXX"
   save_CC="$CC"
   save_LIBS="$LIBS"
   compiler_subst () { subst="[$]1"; shift; orig="[$]1"; shift; echo "[$]subst [$][@]"; };

   AC_LANG_PUSH([C])
   AX_MPI([have_mpicc="yes"],
          [have_mpicc="no";
           AC_MSG_WARN([No C MPI implementation found: MPI example C programs will not be built.])])
   AS_IF([test AS_VAR_GET([have_mpicc]) != no],[
     CC=`compiler_subst $MPICC $CC`
     LIBS="$MPILIBS $LIBS"
   ],[
     CC="$save_CC"
     LIBS="$save_LIBS"
   ])
   AC_SUBST(MPICC)
   AC_LANG_POP([C])

   AC_LANG_PUSH([C++])
   AX_MPI([have_mpicxx="yes"],
          [have_mpicxx="no";
           AC_MSG_WARN([No C++ MPI implementation found: MPI example C++ programs will not be built.])])
   AS_IF([test AS_VAR_GET([have_mpicxx]) != no],[
     CXX=`compiler_subst $MPICXX $CXX`
     LIBS="$MPILIBS $LIBS"
   ],[
     CXX="$save_CXX"
     LIBS="$save_LIBS"
   ])
   AC_SUBST(MPICXX)
   AC_LANG_POP([C++])
])

# check for needed Boost libraries
BOOST_REQUIRE
BOOST_STATIC
BOOST_SERIALIZATION
BOOST_MPI([], 
          [], 
          [# disable MPI/C++ altogether
           have_mpicxx=no;
           AC_MSG_WARN([Could not find a usable Boost.MPI library; not compiling MPI example C++ programs.])])
BOOST_TEST([], 
           [have_boost_test=yes], 
           [have_boost_test=no;
            AC_MSG_WARN([Could not find a usable Boost.Test library; not compiling test programs.])])
AM_CONDITIONAL([HAVE_BOOST_TEST], [test "x$have_boost_test" != no])


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
