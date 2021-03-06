/**
 * @file   test_rank_int.cpp
 *
 * Correctness test for the ``Rheinfall`` C++ implementation, with `long
 * long` coefficients.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2010, 2011 riccardo.murri@gmail.com.  All rights reserved.
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

// ensure assertions are enabled in test cases
#ifdef NDEBUG
# undef NDEBUG
#endif


#include <config.h>
// ignore GMP even if defined, since we're not going to test it
#ifdef HAVE_GMPXX
# undef HAVE_GMPXX
#endif

typedef long long val_t;
typedef long coord_t;

#include <rank.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>


#define BOOST_TEST_MODULE rheinfall
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE( test_rank_int );
#include "test_rank.cpp.inc"
BOOST_AUTO_TEST_SUITE_END()
