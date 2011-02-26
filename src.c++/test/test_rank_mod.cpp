/**
 * @file   test_rank_long.cpp
 *
 * Correctness test for the `Rheinfall` C++ implementation, with `long
 * double` coefficients.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010 riccardo.murri@gmail.com.  All rights reserved.
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
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

#include <config.h>
// ignore GMP even if defined, since we're not going to test it
#ifdef HAVE_GMPXX
# undef HAVE_GMPXX
#endif

#include <modular.hpp>

typedef long mod_int_t;
typedef modular::Modular<mod_int_t> val_t;
typedef long coord_t;

#include <rheinfall.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>


#define BOOST_TEST_MODULE rheinfall
#include <boost/test/unit_test.hpp>

// test fixture
struct F {
  // setup fixture
  F()  { val_t::global_set_modulus(2038076783); };
  // teardown fixture
  ~F() { };
};


BOOST_FIXTURE_TEST_SUITE( test_rank_mod, F )
#include "test_rank.cpp.inc"
BOOST_AUTO_TEST_SUITE_END()
