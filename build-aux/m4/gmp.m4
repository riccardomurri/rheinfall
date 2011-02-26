# Check for GMP
# Modified by Pascal Giorgi, 2003-12-03
# Modified by Riccardo Murri <riccardo.murri@gmail.com>, 2010-10-16

dnl AX_CHECK_GMP ([MINIMUM-VERSION], [PATH], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl ----------------------------------------------------------------------------------
dnl 
dnl Check if the GNU Multiprecision library is available for compiling and
dnl linking; if successful, define GMP_CPPFLAGS and GMP_LIBS and run shell
dnl commands ACTION-IF_FOUND; otherwise run shell commands
dnl ACTION-IF-NOT-FOUND.  Directories specified in PATH (space-separated
dnl list) are added to the compiler/linker search path while looking for
dnl GMP If MINIMUM-VERSION is given, then require that GMP has at least
dnl that version number; by default, any version is acceptable.
dnl
AC_DEFUN([AX_CHECK_GMP], 
    [
        gmp_req_version=ifelse([$1], ,0.0.0,$1)
        gmp_req_major=`echo $gmp_req_version | cut -d. -f1`
        gmp_req_minor=`echo $gmp_req_version | cut -d. -f2`
        if test "x$gmp_req_minor" = x; then
            gmp_req_minor=0
        fi
        gmp_req_patchlevel=`echo $gmp_req_version | cut -d. -f3`
        if test "x$gmp_req_patchlevel" = x; then
            gmp_req_patchlevel=0
        fi
        
        AC_MSG_CHECKING(for GMP)

        gmp_save_CPPFLAGS="${CPPFLAGS}"
        gmp_save_LDFLAGS="${LDFLAGS}"
        gmp_save_LIBS="${LIBS}"
        gmp_save_LD_LIBRARY_PATH="${LD_LIBRARY_PATH}"

        for GMP_HOME in $2 ${DEFAULT_CHECKING_PATH} ''; do        
            if test -r "$GMP_HOME/include/gmp.h"; then :
                GMP_CPPFLAGS="-I${GMP_HOME}/include"
                GMP_LDFLAGS="-L${GMP_HOME}/lib"
                GMP_LIBS="-lgmp"        
                gmp_LD_LIBRARY_PATH="${GMP_HOME}/lib:"
        
                CPPFLAGS="${gmp_save_CPPFLAGS} ${GMP_CPPFLAGS}"
                LDFLAGS="${gmp_save_LDFLAGS} ${GMP_LDFLAGS}"
                LIBS="${gmp_save_LIBS} ${GMP_LIBS}"
                LD_LIBRARY_PATH="${gmp_LD_LIBRARY_PATH}${LD_LIBRARY_PATH}" # XXX: is this portable??
                export LD_LIBRARY_PATH
        
                AC_LANG_PUSH([C])
                AC_TRY_LINK(
                    [#include <gmp.h>],
                    [mpz_t a; mpz_init (a);],
                    [ # check for header OK
                        AC_MSG_RESULT(found)
                        AC_MSG_CHECKING(whether GMP is version $gmp_req_version or greater)
                        AC_TRY_RUN(
                            [
#include <gmp.h>
int main() {
    if (__GNU_MP_VERSION > ]$gmp_req_major[) return 0;
    else if (__GNU_MP_VERSION == ]$gmp_req_major[)
        if (__GNU_MP_VERSION_MINOR > ]$gmp_req_minor[) return 0;
        else if (__GNU_MP_VERSION_MINOR == ]$gmp_req_minor[)
#ifndef __GNU_MP_VERSION_PATCHLEVEL
#define __GNU_MP_VERSION_PATCHLEVEL 0
#endif
             if (__GNU_MP_VERSION_PATCHLEVEL >= ]$gmp_req_patchlevel[) return 0;
    return 1;
}
                            ],[ # try run version check OK
                                have_gmp="yes"
                                break
                            ],[ # try run version check unsuccessful
                                have_gmp="no"        
                            ],[ # cross-compiler should be catched outside
                                AC_MSG_FAILURE(cross-compilation detected in inner test. Perhaps a bug?)
                            ])
                    ],[ # check for header unsuccessful
                        have_gmp='no'
                    ],[ # configure thinks we're cross compiling
                        have_gmp='yes'
                        AC_MSG_WARN([WARNING: You appear to be cross compiling, so there is no way to determine whether your GMP version is new enough. I am assuming it is.])
                        break
                    ])
                AC_LANG_POP([C])
            else
                have_gmp='no'
            fi
            if test "x$have_gmp" = 'xno'; then :
                gmp_failed_paths="$gmp_failed_paths $GMP_HOME"
                unset GMP_CPPFLAGS
                unset GMP_LIBS
                unset GMP_LDFLAGS
                unset gmp_LD_LIBRARY_PATH
            fi
        done

        AC_MSG_RESULT([$have_gmp])
        if test "x$have_gmp" = "xyes"; then :
            AC_SUBST([GMP_CPPFLAGS], [${GMP_CPPFLAGS}])
            AC_SUBST([GMP_LDFLAGS], [${GMP_LDFLAGS}])
            AC_SUBST([GMP_LIBS], [${GMP_LIBS}])
            AC_DEFINE([HAVE_GMP], 1,
                [Define if GMP C++ bindings are available])
            ifelse([$3], , :, [$3])
        else
            if test -n "$gmp_failed_paths"; then
                AC_MSG_WARN([Could not find GMP >=$gmp_req_version in any of the following directories: "$gmp_failed_paths"])
            fi
            ifelse([$4], , :, [$4])
        fi
        
        CPPFLAGS="${gmp_save_CPPFLAGS}"
        LIBS="${gmp_save_LIBS}"
        LDFLAGS="${gmp_save_LDFLAGS}"
        LD_LIBRARY_PATH="${gmp_save_LD_LIBRARY_PATH}"
    ])



dnl AX_CHECK_GMPXX ([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl ----------------------------------------------------------------------------------
dnl 
dnl Check if the C++ wrapper classes are provided in the GNU Multiprecision library;
dnl if successful, define GMPXX_CPPFLAGS and GMPXX_LIBS and run shell
dnl commands ACTION-IF_FOUND; otherwise run shell commands
dnl ACTION-IF-NOT-FOUND.  Requires that AX_CHECK_GMP has already been run.
dnl
AC_DEFUN([AX_CHECK_GMPXX], 
    [
        AC_MSG_CHECKING(whether GMP was compiled with --enable-cxx)

        gmp_save_CPPFLAGS="${CPPFLAGS}"
        gmp_save_LDFLAGS="${LDFLAGS}"
        gmp_save_LIBS="${LIBS}"
        gmp_save_LD_LIBRARY_PATH="${LD_LIBRARY_PATH}"

        CPPFLAGS="${gmp_save_CXXFLAGS} ${GMP_CPPFLAGS}"
        LIBS="${gmp_save_LIBS} -lgmpxx ${GMP_LIBS}"
        
        # See if GMP was compiled with --enable-cxx
        AC_LANG_PUSH([C++])
        AC_TRY_LINK(
            [#include <gmpxx.h>], [
    mpz_class a(2),b(3),c(5); 
    if (a+b == c) return 0; else return 1; 
            ],[ # try run gmpxx test OK
                have_gmpxx="yes"
            ],[ # try run gmpxx test unsuccessful
                have_gmpxx="no"        
            ],[ # cross-compilation detected
                have_gmpxx='yes'
                AC_MSG_WARN([You appear to be cross compiling, so there is no way to determine whether your GMP has C++ bindings. I am assuming they are available.])
            ])
        AC_LANG_POP([C++])

        AC_MSG_RESULT([$have_gmpxx])
        if test "x$have_gmpxx" = 'yes'; then :
            AC_DEFINE(HAVE_GMPXX,1,
                [Define if GMP C++ bindings are available])
            AC_SUBST([GMPXX_CPPFLAGS], [${GMP_CPPFLAGS}])
            AC_SUBST([GMPXX_LDFLAGS], [${GMP_LDFLAGS}])
            AC_SUBST([GMPXX_LIBS], [-lgmpxx ${GMP_LIBS}])
            ifelse([$1], , :, [$1])
        else
            ifelse([$2], , :, [$2])
        fi
        
        CPPFLAGS="${gmp_save_CPPFLAGS}"
        LDFLAGS="${gmp_save_LDFLAGS}"
        LIBS="${gmp_save_LIBS}"
    ])
