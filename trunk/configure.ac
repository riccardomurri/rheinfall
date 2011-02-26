# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([rheinfall], [0.1], [riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c/rheinfall.c])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
gl_EARLY
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB

# Initialize GNUlib
gl_INIT

# check for MPI
save_CXX="$CXX"
save_CC="$CC"
save_LIBS="$LIBS"
AC_LANG([C])
AX_MPI([HAVE_MPI="yes"],
       [HAVE_MPI="no"
           AC_MSG_WARN([No C MPI implementation found: rheinfall will be built without MPI support.])])
AS_IF([test AS_VAR_GET([HAVE_MPI]) != no],[
  AC_LANG([C++])
  AX_MPI([HAVE_MPI="yes"],
         [HAVE_MPI="no"
             AC_MSG_WARN([No C++ MPI implementation found: rheinfall++/MPI will be built without MPI support.])])
])
compiler_subst () { subst="[$]1"; shift; orig="[$]1"; shift; echo "[$]subst [$][@]"; };
AS_IF([test AS_VAR_GET([HAVE_MPI]) != no],[
  CC=`compiler_subst $MPICC $CC`
  CXX=`compiler_subst $MPICXX $CXX`
  LIBS="$MPILIBS $LIBS"
],[
  CC="$save_CC"
  CXX="$save_CXX"
  LIBS="$save_LIBS"
])
AC_SUBST(HAVE_MPI)

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
        lib/Makefile
])
AC_OUTPUT(Makefile)
