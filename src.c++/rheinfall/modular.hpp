/**
 * @file   modular.hpp
 *
 * Interface of the modular class.
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

#ifndef MODULAR_HPP
#define MODULAR_HPP


#include <iostream>


namespace modular {

  template<typename val_t>
  class Modular 
  {
    val_t val_;
    static val_t modulus_;
  public:
    Modular(const val_t& value) : val_(value) { };
    Modular() : val_(0) { };

    static void global_set_modulus(const val_t& modulus) { modulus_ = modulus; };

    Modular<val_t>& operator=(const Modular<val_t>& other) { val_ = other.val_; return *this; }
    Modular<val_t>& operator=(const int& other) { val_ = other; return *this; }
    Modular<val_t>& operator=(const long& other) { val_ = other; return *this; }
    Modular<val_t>& operator=(const long long& other) { val_ = other; return *this; }

    Modular<val_t>& operator+=(const Modular<val_t>& other) { val_ += other.val_; val_ %= modulus_; return *this; }
    Modular<val_t>& operator-=(const Modular<val_t>& other) { val_ -= other.val_; val_ %= modulus_; return *this; }
    Modular<val_t>& operator*=(const Modular<val_t>& other) { val_ *= other.val_; val_ %= modulus_; return *this; }
    Modular<val_t>& operator/=(const Modular<val_t>& other) { val_ *= other.inverse().val_; val_ %= modulus_; return *this; }


    /// equality testing
    friend bool operator==(const Modular<val_t>& a, const Modular<val_t>& b) { return 0 == (a.val_ - b.val_) % modulus_; }
    friend bool operator==(const Modular<val_t>& a, const val_t& b) { return 0 == (a.val_ - b) % modulus_; }
    friend bool operator==(const val_t& a, const Modular<val_t>& b) { return 0 == (a - b.val_) % modulus_; }
    friend bool operator!=(const Modular<val_t>& a, const Modular<val_t>& b) { return not (a == b); };
    friend bool operator!=(const Modular<val_t>& a, const val_t& b) { return not (a == b); };
    friend bool operator!=(const val_t& a, const Modular<val_t>& b) { return not (a == b); };

    friend Modular<val_t> operator+(const Modular<val_t>& a, const Modular<val_t>& b) { return Modular<val_t>((a.val_ + b.val_) % Modular<val_t>::modulus_); };

    friend Modular<val_t> operator-(const Modular<val_t>& a, const Modular<val_t>& b) { return Modular<val_t>((Modular<val_t>::modulus_ + a.val_ - b.val_) % Modular<val_t>::modulus_); };
    
    friend Modular<val_t> operator*(const Modular<val_t>& a, const Modular<val_t>& b) { return Modular<val_t>((a.val_ * b.val_) % Modular<val_t>::modulus_); };

    friend Modular<val_t> operator/(const Modular<val_t>& a, const Modular<val_t>& b) { return a * b.inverse(); };

    friend std::ostream& operator<<(std::ostream& out, Modular<val_t> a) { out << a.val_ << "[mod " << Modular<val_t>::modulus_ << "]"; return out; };

    /// unary minus: additive inverse
    friend Modular<val_t> operator-(const Modular<val_t>& a) { return Modular<val_t>((Modular<val_t>::modulus_ - a.val_) % Modular<val_t>::modulus_); };
    
    /// multiplicative inverse
    Modular<val_t> inverse() const
    {
      // compute u,v such that u*val_ + v*modulus_ = 1
      register val_t u, v, a, b;
      u = 1; a = val_;
      v = 0; b = modulus_;
      while (b != 0)
        {
          register val_t q, t1, t2;
          q = a / b;
          t1 = u - q * v; 
          t2 = a - q * b;
          u = v; 
          a = b; 
          v = t1; 
          b = t2;
        }
      return Modular<val_t>(a < 0 ? -u : u);
    };
      
  }; // class Modular<...>
    
  template<typename val_t> val_t Modular<val_t>::modulus_;

}; // namespace modular


#endif // MODULAR_HPP
