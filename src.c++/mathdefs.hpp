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

#include <cstdlib> // std::abs


namespace rheinfall {

  //
  // the GCD functions are only needed for integral types
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
  val_t gcd(val_t const n, val_t const m) 
  { 
    return _gcd(std::abs(n), std::abs(m)); 
  };

  /** Given the leading entries of two rows X and Y, set @a a and @a b
      so that a*X+b*Y is a row whose first non-zero entry is strictly
      to the right of the leading entries of X and Y. The generic
      templated definition is correct for the integers and other PIDs. */
  template <typename val_t>
  static inline
  void get_row_multipliers(val_t const lead_x, val_t const lead_y, 
                           val_t& a, val_t& b) 
  {
    assert(lead_x != 0);
    const val_t c = gcd(lead_x, lead_y); // guaranteed >= 0
    a = - lead_y / c;
    b = lead_x / c;
  };

  //
  // override the generic def for `double` and other pseudo-R types
  //
#define SPECIALIZE_GET_ROW_MULITPLIERS_FOR_FIELD(T) \
  template <>                                       \
  void get_row_multipliers<T>(T const lead_x,       \
                              T const lead_y,       \
                              T& a, T& b)           \
  {                                                 \
    assert(lead_x != 0);                            \
    a = - lead_y / lead_x;                          \
    b = 1;                                          \
  };                                                \

  SPECIALIZE_GET_ROW_MULITPLIERS_FOR_FIELD(float)
  SPECIALIZE_GET_ROW_MULITPLIERS_FOR_FIELD(double)
#ifdef HAVE_LONG_DOUBLE
  SPECIALIZE_GET_ROW_MULITPLIERS_FOR_FIELD(long double)
#endif

}; // namespace rheinfall


#endif // RF_MATHDEFS_HPP
