# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([rheinfall], [0.1], [riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c/rheinfall.c])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
gl_EARLY
AC_PROG_CC_C99
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB

# Initialize GNUlib
gl_INIT

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
