/**
 * @file   xint.hpp
 *
 * Declaration and definition of mathematical operations on classes 
 * from the Bost.XInt multiprecision integer library.
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

#ifndef RF_TYPES_XINT_HPP
#define RF_TYPES_XINT_HPP


#include <boost/xint/integer.hpp>

#include <rheinfall/types.hpp>


namespace rheinfall {
  // mpz_t -- Integer ring
  
  /** Adapt Boost.Xint's @c boost::xint::integer to our calling convention. */
  static inline
  void rf_xint_gcd(boost::xint::integer &result, 
                   boost::xint::integer const &p, 
                   boost::xint::integer const &q)
  {
    result = boost::xint::gcd(p, q);
  };
  RF_TYPE_IS_GCD_RING(boost::xint::integer, rf_xint_gcd)
  
  // use in-place update in Row::gaussian_elimination()
  template<> class use_inplace_update<boost::xint::integer> : public mpl::bool_<true> { };

}; // namespace rheinfall


#endif // RF_TYPES_XINT_HPP
