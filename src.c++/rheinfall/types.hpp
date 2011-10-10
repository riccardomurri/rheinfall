/**
 * @file   types.hpp
 *
 * Declaration and definition of mathematical operations on supported
 * types.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * See http://www.gnu.org/licenses/gpl.html for licence details.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */

#ifndef RF_TYPES_HPP
#define RF_TYPES_HPP


#include <boost/mpl/bool.hpp>
namespace mpl = boost::mpl;
#include <boost/static_assert.hpp>


#include <algorithm>
#include <cassert>
#include <cstdlib> // std::abs
#include <iostream>
#include <limits>
#include <sstream>


namespace rheinfall {

  /** Return @c true if term @a lead_x is a better pivot for Gaussian
      Elimination than @a lead_y. Default implementation returns @c
      false so that rows are not exchanged; see `Processor::step()`. */
  template <typename val_t>
  static inline
  bool first_is_better_pivot(val_t const& lead_x, val_t const& lead_y)
  {
    return false;
  }


  /** Return @c true if value @a candidate satisfies the "threshold
      pivoting" constraint, given the "best" pivot value @a best and
      the threshold @a ths.  The default implementation always returns
      @c false so that any pivot is considered a candidate; see
      `Processor::step()`.  (The reason this is implemented as a
      function here is that the "threshold pivoting" relation is
      different on the integers and on floating-point numbers.) */
  template <typename val_t>
  static inline
  bool good_enough_pivot(val_t const& best, val_t const& ths, 
                         val_t const& candidate)
  {
    return true;
  }


  /** Given the leading entries of two rows X and Y, set @a a and @a b
      so that a*X+b*Y is a row whose first non-zero entry is strictly
      to the right of the leading entries of X and Y. The generic
      templated definition will error out at compile time: you have to
      provide a specialization for any type you wish to use. See
      macros RF_TYPE_IS_DIVISION_RING and RF_TYPE_IS_GCD_RING */
  template <typename val_t>
  static inline
  void get_row_multipliers(val_t const& lead_x, val_t const& lead_y,
                           val_t& a, val_t& b)
  {
    // see http://www.boost.org/doc/libs/1_43_0/doc/html/boost_staticassert.html
    BOOST_STATIC_ASSERT(sizeof(val_t) == 0);
  };



  /** @def RF_TYPE_IS_GCD_RING
   *
   * Declare template specializations for instances of a type @c T
   * modeling a ring with a GCD (greatest common divisor) function.
   * Second argument @c GCDFUNC is a function @c {void gcd(T& a, T& p, T& q)}
   * that computes the GCD of @c p @c q and stores it into @c a;
   * it must also guarantee that @c a >= 0.
   */
#define RF_TYPE_IS_GCD_RING(T, GCDFUNC)             \
  template <>                                       \
  void get_row_multipliers< T >(T const& lead_x,    \
                                T const& lead_y,    \
                                T& a, T& b)         \
  {                                                 \
    assert(lead_x != 0);                            \
    T c;                                            \
    GCDFUNC(c, lead_x, lead_y);                     \
    assert(c >= 0);                                 \
    a = - lead_y / c;                               \
    b = lead_x / c;                                 \
  };                                                \
                                                    \
  template <>                                       \
  bool first_is_better_pivot<T>(T const& lead_x,    \
                                T const& lead_y)    \
  {                                                 \
     /* Keep absolute value low,                    \
      * so GCD is quicker to find.                  \
      * This implementation is convoluted so to     \
      * 1. Avoid calling `abs` on `mpz_t` & others  \
      * 2. In POD types, we can always take the     \
      *    negative of a positive number, but the   \
      *    converse is not necessarily true.        \
      */                                            \
  return                                            \
    ((lead_x > 0) ?                                 \
     ((lead_y > 0) ?                                \
      (lead_x < lead_y)     /* x,y > 0  */          \
      : (-lead_x > lead_y)) /* x>0, y<0 */          \
     :                                              \
     ((lead_y > 0) ?                                \
      (lead_x > -lead_y)    /* x<0, y>0 */          \
      : (lead_x > lead_y)));  /* x,y < 0  */        \
  };                                                \
                                                    \
  template <>                                       \
  bool good_enough_pivot(T const& best,             \
                         T const& ths,              \
                         T const& candidate)        \
  {                                                 \
    return (candidate <= (best*ths));               \
  };                                                \

  /// The GCD function for integral types.
  //

  /** GCD implementation for basic integral types. */
  template <typename val_t>
  static inline
  val_t _gcd(val_t const n, val_t const m)
  {
    return m==0? n : _gcd(m, n%m);
  };

  /** Return GCD of the given arguments.  Result is always a non-negative integer. */
  template <typename val_t>
  static inline
  void pod_gcd(val_t &result, val_t const &n, val_t const &m)
  {
    result = _gcd<val_t>(std::abs(n), std::abs(m));
  };

  RF_TYPE_IS_GCD_RING(int,       pod_gcd<int>)
  RF_TYPE_IS_GCD_RING(long,      pod_gcd<long>)
#ifdef HAVE_LONG_LONG_INT
  RF_TYPE_IS_GCD_RING(long long, pod_gcd<long long>)
#endif



  /** @def RF_TYPE_IS_DIVISION_RING
   * 
   * Declare template specializations for instances of a type
   * modeling a division ring (i.e., all four algebraic operations
   * supported, @c operator/ is inverse to @c operator* but @c
   * operator* is not necessarily commutative).
   */
#define RF_TYPE_IS_DIVISION_RING(T)                 \
  template <>                                       \
  void get_row_multipliers<T>(T const& lead_x,      \
                              T const& lead_y,      \
                              T& a, T& b)           \
  {                                                 \
    assert(lead_x != 0);                            \
    a = - lead_y / lead_x;                          \
    b = 1;                                          \
  };                                                \
                                                    \
  template <>                                       \
  bool first_is_better_pivot<T>(T const& lead_x,    \
                                T const& lead_y)    \
  {                                                 \
     /* Keep absolute value as high as possible,    \
      * to reduce numeric error.                    \
      * This implementation is convoluted so to     \
      * avoid calling `abs` on `mpz_f` & others.    \
      */                                            \
  return                                            \
    ((lead_x > 0) ?                                 \
     ((lead_y > 0) ?                                \
      (lead_x > lead_y)     /* x,y > 0  */          \
      : (-lead_x < lead_y)) /* x>0, y<0 */          \
     :                                              \
     ((lead_y > 0) ?                                \
      (lead_x < -lead_y)    /* x<0, y>0 */          \
      : (lead_x < lead_y)));  /* x,y < 0  */        \
  };                                                \
                                                    \
  template <>                                       \
  bool good_enough_pivot(T const& best,             \
                         T const& ths,              \
                         T const& candidate)        \
  {                                                 \
    return (candidate >= (best*ths));               \
  };                                                \

  RF_TYPE_IS_DIVISION_RING(float)
  RF_TYPE_IS_DIVISION_RING(double)
#ifdef HAVE_LONG_DOUBLE
  RF_TYPE_IS_DIVISION_RING(long double)
#endif



  /** @def RF_TYPE_IS_UNORDERED_DIVISION_RING
   * 
   * Declare template specializations for instances of a type
   * modeling an unordered division ring: i.e., like @ref
   * RF_TYPE_IS_DIVISION_RING but no relational operators are
   * defined, so fall back to the generic implementation of @ref
   * first_is_better_pivot<val_t>().
   */
#define RF_TYPE_IS_UNORDERED_DIVISION_RING(T)       \
  template <>                                       \
  void get_row_multipliers< T >(T const& lead_x,    \
                                T const& lead_y,    \
                                T& a, T& b)         \
  {                                                 \
    assert(lead_x != 0);                            \
    a = - lead_y / lead_x;                          \
    b = 1;                                          \
  };                                                \



  /** @def RF_TYPE_IS_MODULAR_RING
   * 
   * Declare template specializations for instances of a type
   * modeling a modular ring, i.e., take advantage of the fact that
   * multiplication is cheaper than division and never overflows.
   */
#define RF_TYPE_IS_MODULAR_RING(T)                  \
  template <>                                       \
  void get_row_multipliers< T >(T const& lead_x,    \
                                T const& lead_y,    \
                                T& a, T& b)         \
  {                                                 \
    assert(lead_x != 0);                            \
    a = - lead_y;                                   \
    b = lead_x;                                     \
  };                                                \


  /** Tag class to determine whether to update DenseRow storage
      in-place in the DenseRow::gaussian_elimination() methods. By
      default, Rheinfall does @a not use in-place update. Specialize
      to an instance of @c boost::mpl::true_ to enable in-place update. */
  template<typename T> class use_inplace_update : public mpl::bool_<false> { };



  /** Compare a value to zero; by default this is an exact comparison. */
  template <typename T> 
  bool is_zero (const T& value) { return (0 == value); };

  /** Traits class for defining the tolerance for zero equality of an
      inexact type. */
  template <typename T>
  struct is_zero_traits 
  {
    /** Maximum positive value that will be be considered equal to zero. */
    static T tolerance;
  };
  template<typename T> T is_zero_traits<T>::tolerance;

  /** @def INEXACT_IS_ZERO
   *
   *  Specialize the @c is_zero template for inexact arithmetic
   *  types. Consider an inexact value to be zero iff its absolute
   *  value is within a certain distance from the exact zero. 
   */
#define INEXACT_IS_ZERO(T)                                     \
  template<>                                                   \
  bool is_zero<T>(const T& value) {                            \
    return (value > -(is_zero_traits<T>::tolerance))           \
            and (value < is_zero_traits<T>::tolerance);        \
  };                                                           \

  INEXACT_IS_ZERO(float)
  INEXACT_IS_ZERO(double)
#ifdef HAVE_LONG_DOUBLE
  INEXACT_IS_ZERO(long double)
#endif

}; // namespace rheinfall


#endif // RF_TYPES_HPP
