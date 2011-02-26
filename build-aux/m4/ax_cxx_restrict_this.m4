# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_cxx_restrict_this.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_RESTRICT_THIS
#
# DESCRIPTION
#
#   Determine whether the C++ compiler supports defining a member
#   function with a restricted "this" pointer.  Define "restrict_this"
#   to the correct spelling; use like this:
#      
#     T::fn() restrict_this { /* code */ }
#
#   Otherwise, define "restrict_this" to be empty.
#
#   Note: the syntax above is a GCC extension. If your C++ compiler
#   has a different way of applying the 'restricted' qualifier to the
#   "this" pointer, please consider reporting it.
#
# LICENSE
#
#   Copyright (c) 2010 Riccardo Murri <riccardo.murri@gmail.com>
#
#   This file includes code originally found in the AC_C_RESTRICT
#   macro definition from the autoconf-2.68 sources.  Copying and
#   distribution of this file are permitted according to the same
#   terms as the autoconf-2.68 source files.  This file is offered
#   as-is, without any warranty.

#serial 1

AC_DEFUN([AX_CXX_RESTRICT_THIS],
[AC_CACHE_CHECK([whether C++ supports GCC's restrict "this" syntax], ax_cv_cxx_restrict_this,
  [ax_cv_cxx_restrict_this=no
   AC_LANG_PUSH([C++])
   # The order here caters to the fact that C++ does not require restrict.
   for ac_kw in __restrict __restrict__ _Restrict restrict; do
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [[#define restrict_this $ac_kw
class Foo {
  private:
      int x_;    
  public:
      Foo(int x) : x_(x) { };
      int bar(Foo* other);
};
int Foo::bar(Foo* other) restrict_this { return this->x_ + other->x_; }
      ]],
      [[
Foo a(3), b(5);
return (8 == a.bar(&b));      
      ]])],
      [ax_cv_cxx_restrict_this=$ac_kw])
     test "$ax_cv_cxx_restrict_this" != no && break
   done
   AC_LANG_POP([C++])
  ])
 AH_VERBATIM([restrict_this],
[/* Define to the keyword(s) used to specify that a member function's
    "this" pointer is unaliased.  Define to nothing if this is not
    supported. */
#undef restrict_this])
 case $ax_cv_cxx_restrict_this in
   no) AC_DEFINE([restrict_this], []) ;;
   *)  AC_DEFINE_UNQUOTED([restrict_this], [$ax_cv_cxx_restrict_this]) ;;
 esac
])# AC_C_RESTRICT


