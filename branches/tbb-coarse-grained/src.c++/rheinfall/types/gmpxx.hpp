/**
 * @file   gmpxx.hpp
 *
 * Declaration and definition of mathematical operations on classes 
 * from the GNU MP library.
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

#ifndef RF_TYPES_GMPXX_HPP
#define RF_TYPES_GMPXX_HPP


# include <gmpxx.h>

#include <rheinfall/types.hpp>

namespace rheinfall {
  // mpz_t -- Integer ring
  
  /** Adapt GMPXX's @c mpz_gcd function to our calling convention. */
  static inline
  void rf_mpz_gcd(mpz_class &result, mpz_class const &p, mpz_class const &q)
  {
    mpz_gcd(result.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());
  };
  RF_TYPE_IS_GCD_RING(mpz_class, rf_mpz_gcd)
  
  // use in-place update in Row::gaussian_elimination()
  template<> class use_inplace_update<mpz_class> : public mpl::bool_<true> { };
  
  
  // mpq_t -- Rational field (representation as integer quotient)
  
  RF_TYPE_IS_DIVISION_RING(mpq_class)
  
  template<> class use_inplace_update<mpq_class> : public mpl::bool_<true> { };
  
  
  // mpf_t -- Rational field (representation as infinite-precision decimal fractional number)
  RF_TYPE_IS_DIVISION_RING(mpf_class)
  
  template<> class use_inplace_update<mpf_class> : public mpl::bool_<true> { };
  
}; // namespace rheinfall

  
#endif // RF_TYPES_GMPXX_HPP
