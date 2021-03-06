/**
 * @file   gmpxx_boost_serialization.hpp
 *
 * Support for serialization of GMP data types with the
 * Boost.Serialization library.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2012 riccardo.murri@gmail.com. All rights reserved.
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

#ifndef GMPXX_BOOST_SERIALIZATION_HPP
#define GMPXX_BOOST_SERIALIZATION_HPP

#include <gmpxx.h>

// Boost.Serialization support for GMP classes
# include <boost/serialization/split_free.hpp>

BOOST_SERIALIZATION_SPLIT_FREE(__mpz_struct)
BOOST_SERIALIZATION_SPLIT_FREE(__mpf_struct)

namespace boost {
namespace serialization {

  // mpz_class
  template<typename Archive>
  inline
  void serialize(Archive& ar, mpz_class& z, const boost::serialization::version_type& version)
  {
    ar & *(z.get_mpz_t()); // `.mp` is the sole data member of mp*_class types
  };

  template<typename Archive>
  inline
  void load(Archive& ar, __mpz_struct& z, const boost::serialization::version_type& version)
  {
    // save original size; we need it to call mp_free_fn
    mp_size_t original_limb_no = z._mp_alloc;
    // do the POD part
    ar & z._mp_alloc & z._mp_size;
    // get ptrs to GMP memory management functions
    void* (*alloc_func_ptr) (size_t);
    void  (*free_func_ptr) (void *, size_t);
    mp_get_memory_functions(&alloc_func_ptr, NULL, &free_func_ptr);
    // free old memory
    if (original_limb_no != 0)
      (*free_func_ptr)(z._mp_d, original_limb_no);
    // make room for arriving limbs
    z._mp_d = static_cast<mp_limb_t*>((*alloc_func_ptr)(z._mp_alloc));
    // now serialize the limbs array
    mp_limb_t *limbp = z._mp_d;
    for (int i = 0; i < z._mp_alloc; ++i) {
      ar & (*limbp);
      ++limbp;
    };
  };

  template<typename Archive>
  inline
  void save(Archive& ar, const __mpz_struct& z,
            const boost::serialization::version_type& version)
  {
    ar & z._mp_alloc & z._mp_size;
    // now serialize the limbs array
    mp_limb_t *limbp = z._mp_d;
    for (int i = 0; i < z._mp_alloc; ++i) {
      ar & (*limbp);
      ++limbp;
    };
  };


  // mpq_class
  template<typename Archive>
  inline
  void serialize(Archive& ar, mpq_class& q, const boost::serialization::version_type& version)
  {
    ar & *(q.get_mpq_t()); // `.mp` is the sole data member of mp*_class types
  };

  template<typename Archive>
  inline
  void serialize(Archive& ar, __mpq_struct& q, const boost::serialization::version_type& version)
  {
    ar & q._mp_num & q._mp_den;
  };


  // mpf_class
  template<typename Archive>
  inline
  void serialize(Archive& ar, mpf_class& f, const boost::serialization::version_type& version)
  {
    ar & *(f.get_mpf_t()); // `.mp` is the sole data member of mp*_class types
  };

  template<typename Archive>
  inline
  void load(Archive& ar, __mpf_struct& f, const boost::serialization::version_type& version)
  {
    // save original size; we need it to call mp_free_fn
    mp_size_t original_limb_no = 1 + f._mp_prec;
    // do the POD part
    ar & f._mp_prec & f._mp_size & f._mp_exp;
    // get ptrs to GMP memory management functions
    void *(*alloc_func_ptr) (size_t);
    void *(*free_func_ptr) (void *, size_t);
    mp_get_memory_functions(&alloc_func_ptr, NULL, &free_func_ptr);
    // free old memory
    if (original_limb_no != 0)
      (*free_func_ptr)(f._mp_d, original_limb_no);
    // make room for arriving limbs
    f._mp_d = static_cast<mp_limb_t*>((*alloc_func_ptr)(1 + f._mp_prec));
    // now serialize the limbs array
    mp_limb_t *limbp = f._mp_d;
    for (int i = 0; i < 1 + f._mp_prec; ++i) {
      ar & (*limbp);
      ++limbp;
    };
  };

  template<typename Archive>
  inline
  void save(Archive& ar, const __mpf_struct& f,
            const boost::serialization::version_type& version)
  {
    ar & f._mp_prec & f._mp_size & f._mp_exp;
    // now serialize the limbs array
    mp_limb_t *limbp = f._mp_d;
    for (int i = 0; i < 1 + f._mp_prec; ++i) {
      ar & (*limbp);
      ++limbp;
    };
  };

} // namespace serialization
} // namespace boost



#endif // GMPXX_BOOST_SERIALIZATION_HPP
