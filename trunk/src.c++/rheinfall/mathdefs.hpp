/**
 * @file   mathdefs.hpp
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

#ifndef RF_MATHDEFS_HPP
#define RF_MATHDEFS_HPP


#include <boost/static_assert.hpp>

#ifdef HAVE_GMPXX
# include <gmpxx.h>
#endif

#include <cassert>
#include <cstdlib> // std::abs



namespace rheinfall {

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
#define RF_TYPE_IS_GCD_RING(T, gcd)                 \
  template <>                                       \
  void get_row_multipliers(T const& lead_x,         \
                           T const& lead_y,         \
                           T& a, T& b)              \
  {                                                 \
    assert(lead_x != 0);                            \
    T c;                                            \
    gcd(c, lead_x, lead_y);                         \
    assert(c >= 0);                                 \
    a = - lead_y / c;                               \
    b = lead_x / c;                                 \
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
  RF_TYPE_IS_GCD_RING(long long, pod_gcd<long long>)

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
#define RF_TYPE_IS_DIVISION_RING(T) \
  template <>                                       \
  void get_row_multipliers<T>(T const& lead_x,      \
                              T const& lead_y,      \
                              T& a, T& b)           \
  {                                                 \
    assert(lead_x != 0);                            \
    a = - lead_y / lead_x;                          \
    b = 1;                                          \
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

}; // namespace rheinfall


#endif // RF_MATHDEFS_HPP
