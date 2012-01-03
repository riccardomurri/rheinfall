## Process this file to produce a configure script:
##   aclocal -I build-aux/m4 && autoheader && autoconf && automake
##

AC_PREREQ([2.65])
AC_INIT([rheinfall],[0.15.10],[riccardo.murri@gmail.com])
AC_CONFIG_SRCDIR([src.c++/rheinfall/rank.hpp])
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


# Check for Intel TBB
AC_ARG_WITH([tbb],
    [AS_HELP_STRING([--with-tbb=yes|no|@<:@path@:>@],
            [Build example programs requiring the Intel TBB library. 
                If argument is 'yes' or empty, the library must be 
                reachable with the standard compiler search path; 
                otherwise the TBB include files have to be in 
                <path>/include and the compiled library has to 
                reside in <path>/lib.  If the argument is 'no', then
                example programs that make use of TBB will not be built.])],
    [ # option '--with-tbb' given
      if test "_$withval" != _no ; then :
          AC_LANG_PUSH([C++])
         save_CPPFLAGS="$CPPFLAGS"
         save_CXXFLAGS="$CXXFLAGS"
         save_LDFLAGS="$LDFLAGS"
         save_LIBS="$LIBS"
         if test "_$withval" != "_yes" ; then :
             TBB_CPPFLAGS="-I$withval/include"
         fi
         CPPFLAGS="$save_CPPFLAGS $TBB_CPPFLAGS"
         AC_CHECK_HEADER([tbb/tbb.h], 
           [ # 'tbb/tbb.h' found
             AC_MSG_CHECKING([for the Intel TBB library])
             # is there a TBB compiled specifically for this arch/compiler?
             tbb_root="$withval"
             if test -r "$tbb_root"/Makefile; then :
                 tbb_arch=$(cd "$tbb_root"; make info | grep '^arch=' | cut -d= -f2)
                 tbb_runtime=$(cd "$tbb_root"; make info | grep '^runtime=' | cut -d= -f2)
                 tbb_libdir_guess="$tbb_root/lib/$tbb_arch/$tbb_runtime"
             fi
             # directory naming scheme in Intel's 'tbb*_lin.tgz' distribution is
             # "lib/$tbb_arch/$tbb_runtime", where $tbb_runtime has the the form
             # "cc<version>_libc<version>_kernel<version>"
             tbb_libdir_other=$(find "$withval/lib/$tbb_arch" -type d -name 'cc*')
             # the '' directory is special-cased below: it looks for
             # TBB in the default compiler search path
             for libdir in $tbb_libdir_guess $tbb_libdir_other ''; do
               TBB_LIBS="-ltbb"
               if test "_$libdir" != '_' ; then :
                 TBB_LDFLAGS="-L$libdir"
               else   
                 TBB_LDFLAGS=""
               fi
               LDFLAGS="$save_LDFLAGS $TBB_LDFLAGS"
               LIBS="$save_LIBS $TBB_LIBS"
               AC_LINK_IFELSE(
                   [AC_LANG_SOURCE([
#include "tbb/task_scheduler_init.h"
using namespace tbb;
int main() {
  task_scheduler_init init;
  return 0;
}
                                ])],
                   [tbb_found=yes;],
                   [tbb_found=no;])
               if test "_$tbb_found" = '_yes' ; then :
                 break
               fi
             done;
             if test "_$tbb_found" = '_yes' ; then :
               AC_MSG_RESULT([TBB library found in "$libdir"])
               have_tbb=yes
               AC_SUBST([HAVE_TBB],     [1])
               AC_SUBST([TBB_CPPFLAGS], ["$TBB_CPPFLAGS"])
               AC_SUBST([TBB_CXXFLAGS], ["-DRF_USE_TBB"])dnl XXX: Rheinfall-specific
               AC_SUBST([TBB_LDFLAGS],  ["$TBB_LDFLAGS"])
               AC_SUBST([TBB_LIBS],     ["$TBB_LIBS"])
             else
               AC_MSG_RESULT([no TBB library found])
               have_tbb=no
               AC_SUBST([TBB_CXXFLAGS], [])
               AC_SUBST([TBB_LDFLAGS],  [])
               AC_SUBST([TBB_LIBS],     [])
               AC_MSG_WARN([Cannot link with the TBB library: not compiling example TBB programs.])
             fi
           ], 
           [ # no include files, disable TBB altogether
             AC_MSG_RESULT([no]) 
             have_tbb=no;        
             AC_MSG_WARN([Cannot include TBB library headers: not compiling example TBB programs.])
           ])
         AC_LANG_POP([C++])
         CPPFLAGS="$save_CPPFLAGS"
         CXXFLAGS="$save_CXXFLAGS"
         LDFLAGS="$save_LDFLAGS"
         LIBS="$save_LIBS"
      fi # $withval != 'no
    ])
AM_CONDITIONAL([HAVE_TBB], [test "_$have_tbb" = '_yes'])


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
