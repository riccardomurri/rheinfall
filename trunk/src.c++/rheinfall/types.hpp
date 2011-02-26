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
 * Copyright (c) 2010 riccardo.murri@gmail.com. All rights reserved.
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

#ifdef HAVE_GMPXX
# include <gmpxx.h>
#endif

#include <cassert>
#include <cstdlib> // std::abs
#include <limits>



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


  /** Given the leading entries of two rows X and Y, set @a a and @a b
      so that a*X+b*Y is a row whose first non-zero entry is strictly
      to the right of the leading entries of X and Y. The generic
      templated definition errors out, you have to provide a
      specialization for any type you wish to use. See macros
      RF_TYPE_IS_DIVISION_RING and RF_TYPE_IS_GCD_RING */
  template <typename val_t>
  static inline
  void get_row_multipliers(val_t const& lead_x, val_t const& lead_y,
                           val_t& a, val_t& b)
  {
    // see http://www.boost.org/doc/libs/1_43_0/doc/html/boost_staticassert.html
    BOOST_STATIC_ASSERT(sizeof(val_t) == 0);
  };


  /** Declare template specializations for instances of a type @c T
   * modeling a ring with a GCD (greatest common divisor) function.
   * Second argument @c gcd is a function @c{void gcd(T& a, T& p, T& q)}
   * that computes the GCD of @c p @c q and stores it into @c{a};
   * it must also guarantee that @c{a >= 0}..
   */
#define RF_TYPE_IS_GCD_RING(T, GCD)                 \
  template <>                                       \
  void get_row_multipliers<T>(T const& lead_x,      \
                              T const& lead_y,      \
                              T& a, T& b)           \
  {                                                 \
    assert(lead_x != 0);                            \
    T c;                                            \
    GCD(c, lead_x, lead_y);                         \
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

  // the GCD function for POD integral types
  //
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

#ifdef HAVE_GMPXX
  /** Adapt GMPXX's @c mpz_gcd function to our calling convention. */
  static inline
  void rf_mpz_gcd(mpz_class &result, mpz_class const &p, mpz_class const &q)
  {
    mpz_gcd(result.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());
  };
  RF_TYPE_IS_GCD_RING(mpz_class, rf_mpz_gcd)
#endif 


  /** Declare template specializations for instances of a type
   *  modeling a division ring (i.e., all four algebraic operations
   *  supported, @c operator/ is inverse to @operator* but @c
   *  operator* is not necessarily commutative).
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

  RF_TYPE_IS_DIVISION_RING(float)
  RF_TYPE_IS_DIVISION_RING(double)
#ifdef HAVE_LONG_DOUBLE
  RF_TYPE_IS_DIVISION_RING(long double)
#endif

#ifdef HAVE_GMPXX
  RF_TYPE_IS_DIVISION_RING(mpq_class)
  RF_TYPE_IS_DIVISION_RING(mpf_class)
#endif 



  /** Implementation of the @a is_zero(const val_t&) function. */
  namespace impl {
    // for exact number types, we can just compare with 0
    template<typename val_t> bool is_zero(const val_t& x, mpl::true_,   mpl::true_)  { return static_cast<val_t>(0) == x; }
    template<typename val_t> bool is_zero(const val_t& x, mpl::true_,  mpl::false_)  { return (x == 0); }
    
    // for std floating point types, we can check if number is within a certain tolerance
    // XXX: it should be possible to choose the tolerance!
    template<typename val_t> bool is_zero(const val_t& x, mpl::false_,  mpl::true_)  
    { 
      const val_t b = std::numeric_limits<val_t>::epsilon() * static_cast<val_t>(2);
      return (-b < x) and (x < b);
    }
    
    // no default implementation for non-std inexact types
    template<typename val_t> bool is_zero(const val_t& x, mpl::false_, mpl::false_)  { BOOST_STATIC_ASSERT(sizeof(val_t) == 0); }
  };


  /** Return @c true if @a x is to be considered zero.  
   * The default implementation does the following:
   * - for @a exact numeric types, it does a direct comparison with 0;
   * - for inexact standard types (i.e., primitive floating-point types), 
   *   @c true is returned if the value @a x falls within range -epsilon .. +epsilon;
   * - for other inexact types, a compile-time assertion is thrown, 
   *   so you should specialize this template with the appropriate function 
   *   for user-defined inexact types.
   * 
   */
  template <typename val_t>
  static inline
  bool is_zero(const val_t& x)
  {
    return impl::is_zero(x, 
               mpl::bool_< 
                   std::numeric_limits<val_t>::is_exact 
               >(),
               mpl::bool_<
                   std::numeric_limits<val_t>::is_specialized 
               >());
  };

#define RF_USE_EXACT_IS_ZERO(T) \
  template<> bool is_zero<T>(const T& x) { return x == 0; }

#ifdef HAVE_GMPXX
  RF_USE_EXACT_IS_ZERO(mpf_class)
  RF_USE_EXACT_IS_ZERO(mpq_class)
  RF_USE_EXACT_IS_ZERO(mpz_class)
#endif


  // by default, use in-place update within DenseRow::gaussian_elimination
  template<typename T> class use_inplace_update : public mpl::bool_<false> { };

  // don't use in-place update for POD types, might give some
  // advantage due to the split of read- and write- accesses to memory
  template<> class use_inplace_update<int> : public mpl::bool_<true> { };
  template<> class use_inplace_update<long> : public mpl::bool_<true> { };
#ifdef HAVE_LONG_LONG_INT
  template<> class use_inplace_update<long long> : public mpl::bool_<true> { };
#endif
  template<> class use_inplace_update<float> : public mpl::bool_<true> { };
  template<> class use_inplace_update<double> : public mpl::bool_<true> { };
#ifdef HAVE_LONG_DOUBLE
  template<> class use_inplace_update<long double> : public mpl::bool_<true> { };
#endif
  
}; // namespace rheinfall


#endif // RF_TYPES_HPP