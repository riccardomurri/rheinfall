/**
 * @file   test_sparserow.cpp
 *
 * Unit tests for the `Rheinfall::SparseRow` class.
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

#define BOOST_TEST_MODULE rheinfall::SparseRow
#include <boost/test/included/unit_test.hpp>

#include <sparserow.hpp>
using namespace rheinfall;

#include <map>


BOOST_AUTO_TEST_CASE( check_sparserow_size_null1 )
{
  SparseRow<int,int> s(0, 0, 1);
  BOOST_CHECK_EQUAL( s.size(), 0 );
}


BOOST_AUTO_TEST_CASE( check_sparserow_size_null2 )
{
  SparseRow<int,int> s(0, 100, 1);
  BOOST_CHECK_EQUAL( s.size(), 0 );
}


BOOST_AUTO_TEST_CASE( check_sparserow_size )
{
  SparseRow<int,int> s(0, 100, 1);
  s.set(77,7);
  BOOST_CHECK_EQUAL( s.size(), 1 );
}


BOOST_AUTO_TEST_CASE( check_sparserow_get_set )
{
  SparseRow<int,int> s(0, 100, 1);

  for(int i = 0; i < 10; ++i)
    s.set(i*i, i+1);
  BOOST_CHECK_EQUAL( s.first_nonzero_column(), 0 );

  int j = 0;
  for(int i = 0; i < s.size(); ++i)
    if (j*j == i) {
      BOOST_CHECK_EQUAL(s.get(i), j+1);
      j++;
    }
    else 
      BOOST_CHECK_EQUAL(s.get(i), 0);
}


BOOST_AUTO_TEST_CASE( check_sparserow_ctor_from_iter )
{
  std::map< int,int > entries;
  for(int i = 0; i < 10; ++i)
    entries[i*i] = i+1;

  SparseRow< int,int > *s = 
    SparseRow< int,int >::new_from_range(entries.begin(), entries.end(), 100);
  BOOST_CHECK_EQUAL( s->first_nonzero_column(), 0 );

  int j = 0;
  for(int i = 0; i < s->size(); ++i)
    if (j*j == i) {
      BOOST_CHECK_EQUAL(s->get(i), j+1);
      j++;
    }
    else 
      BOOST_CHECK_EQUAL(s->get(i), 0);
}


BOOST_AUTO_TEST_CASE( check_sparserow_fill_in )
{
  SparseRow<int,int> s(0, 99, 1);
  for (int i = 0; i < 100; ++i) {
    s.set(i, i);
    const double fill_in = s.fill_in();
    BOOST_CHECK( i <= fill_in and fill_in <= i+1 );
  }
}


BOOST_AUTO_TEST_CASE( check_sparserow_adjust )
{
  class S: public SparseRow<int,int> {
  public:
    S(int size) : SparseRow<int,int>::SparseRow(0, size-1, 0) { };
    S* adjust() { return static_cast<S*>(this->SparseRow<int,int>::adjust()); }
  };

  S s(100);
  s.set(77,7);

  S* ss = s.adjust();
  BOOST_CHECK ( ss->first_nonzero_column() == 77 );
  BOOST_CHECK ( ss->size() == 0 );
}


BOOST_AUTO_TEST_CASE( check_sparserow_ge_with_sparse1 )
{
  SparseRow<int,int> *p1 = new SparseRow<int,int>(0,100,2);
  BOOST_CHECK_EQUAL( p1->first_nonzero_column(), 0 );
  BOOST_CHECK_EQUAL( p1->get(0), 2 );

  SparseRow<int,int> *p2 = new SparseRow<int,int>(0,100,5);
  p2->set(77, 7);
  BOOST_CHECK_EQUAL( p2->first_nonzero_column(), 0 );
  BOOST_CHECK_EQUAL( p2->get(0), 5 );
  BOOST_CHECK_EQUAL( p2->get(77), 7 );

  // check elimination
  const int result_size = p2->size() - 1; // p2 will be invalid after GE
  SparseRow<int,int>* p3 = p1->gaussian_elimination(p2);
  BOOST_CHECK ( p1 != p3 );
  BOOST_CHECK ( p2 != p3 );
  BOOST_CHECK_EQUAL( p3->first_nonzero_column(), 77 );
  BOOST_CHECK_EQUAL( p3->size(), result_size );
}


BOOST_AUTO_TEST_CASE( check_sparserow_ge_with_sparse2 )
{
  SparseRow<int,int> *p1 = new SparseRow<int,int>(0,100,2);
  // init s1 - ones at odd places
  for(int i = 1; i <= 100; i += 2)
    p1->set(i, 3);

  SparseRow<int,int> *p2 = new SparseRow<int,int>(0,100,5);
  // init d2 - ones at even places
  for(int i = 2; i <= 100; i += 2)
    p2->set(i, 7);

  // check elimination
  const int result_size = p1->size() + p2->size() - 1; // p2 will be invalid after GE
  SparseRow<int,int>* p3 = p1->gaussian_elimination(p2);
  BOOST_CHECK_EQUAL( p3->size(), result_size );
  BOOST_CHECK_EQUAL( p3->first_nonzero_column(), 1 );
  for (int i = 1; i <= 100; ++i)
    BOOST_CHECK_EQUAL(p3->get(i), (i%2 == 0)? 2*7 : -5*3);
}
