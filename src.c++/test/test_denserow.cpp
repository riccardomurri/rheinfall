
#define BOOST_TEST_MODULE rheinfall::DenseRow
#include <boost/test/included/unit_test.hpp>

#include <denserow.hpp>
using namespace rheinfall;


BOOST_AUTO_TEST_CASE( check_denserow_size_null )
{
  DenseRow<int,int> d(0,0);
  BOOST_CHECK_EQUAL( d.size(), 0 );
}


BOOST_AUTO_TEST_CASE( check_denserow_size )
{
  DenseRow<int,int> d(0,100);
  BOOST_CHECK_EQUAL( d.size(), 100 );
}


BOOST_AUTO_TEST_CASE( check_denserow_get_set )
{
  DenseRow<int,int> d(0,100);

  for(int i = 0; i < 10; ++i)
    d.set(i*i, i+1);
  BOOST_CHECK_EQUAL( d.first_nonzero_column(), 0 );

  int j = 0;
  for(int i = 0; i < d.size(); ++i)
    if (j*j == i) {
      BOOST_CHECK_EQUAL(d.get(i), j+1);
      j++;
    }
    else 
      BOOST_CHECK_EQUAL(d.get(i), 0);
}


BOOST_AUTO_TEST_CASE( check_denserow_adjust )
{
  class D: public DenseRow<int,int> {
  public:
    D(int size) : DenseRow<int,int>::DenseRow(0, size-1) { };
    D* adjust() { return static_cast<D*>(this->DenseRow<int,int>::adjust()); }
  };

  D d(100);
  d.set(77, 7);

  D* dd = d.adjust();
  BOOST_CHECK_EQUAL ( dd->first_nonzero_column(), 77 );
  BOOST_CHECK_EQUAL ( dd->get(77), 7 );
}


BOOST_AUTO_TEST_CASE( check_denserow_ge_with_dense1 )
{
  DenseRow<int,int>*  p1 = new DenseRow<int,int>(0,100);
  p1->set(0, 2);
  BOOST_CHECK_EQUAL( p1->first_nonzero_column(), 0 );

  DenseRow<int,int>*  p2 = new DenseRow<int,int>(0,100);
  p2->set(0, 5);
  p2->set(77, 7);
  BOOST_CHECK_EQUAL( p2->first_nonzero_column(), 0 );
  BOOST_CHECK_EQUAL( p2->get(0), 5 );
  BOOST_CHECK_EQUAL( p2->get(77), 7 );

  // check elimination
  DenseRow<int,int>* p3 = p1->gaussian_elimination(p2);
  BOOST_CHECK ( p1 != p3 );
  BOOST_CHECK ( p2 != p3 );
  BOOST_CHECK_EQUAL( p3->first_nonzero_column(), 77 );
  BOOST_CHECK_EQUAL( p3->get(p3->first_nonzero_column()), 2*7 );
  BOOST_CHECK_EQUAL( p3->size(), p1->size() - 77 );
  for (int i = 0; i < p3->size(); ++i)
    BOOST_CHECK_EQUAL( p3->get(1 + i + p3->first_nonzero_column()), 0 );
}


BOOST_AUTO_TEST_CASE( check_denserow_ge_with_dense2 )
{
  DenseRow<int,int>*  p1 = new DenseRow<int,int>(0,100);
  DenseRow<int,int>*  p2 = new DenseRow<int,int>(0,100);
  
  // init d1 - only at odd places
  p1->set(0, 2);
  for(int i = 1; i <= 100; i += 2)
    p1->set(i, 3);

  // init d2 - only at even places
  p2->set(0, 5);
  for(int i = 2; i <= 100; i += 2)
    p2->set(i, 7);

  // check elimination
  DenseRow<int,int>* p3 = p1->gaussian_elimination(p2);
  BOOST_CHECK_EQUAL( p3->size(), p1->size()-1 );
  BOOST_CHECK_EQUAL( p3->first_nonzero_column(), 1 );
  for (int i = 1; i <= 100; ++i)
    BOOST_CHECK_EQUAL(p3->get(i), (i%2 == 1)? -5*3 : 2*7);
}


//
// interaction with SparseRow
//

#include <sparserow.hpp>

BOOST_AUTO_TEST_CASE( check_denserow_ctor_from_sparserow1 )
{
  SparseRow<int,int> s(0, 100, 1);
  for (int i = 0; i <= 10; ++i)
    s.set(i*i, i+1);

  DenseRow<int,int> d(&s);
  BOOST_CHECK_EQUAL( d.first_nonzero_column(), 0 );
  int j = 0;
  for(int i = 0; i < d.size(); ++i)
    if (j*j == i) {
      BOOST_CHECK_EQUAL(d.get(i), j+1);
      j++;
    }
    else 
      BOOST_CHECK_EQUAL(d.get(i), 0);
}


BOOST_AUTO_TEST_CASE( check_denserow_ctor_from_sparserow2 )
{
  // real example from M05-D5
  SparseRow<int,int> s(0, 7259, 1);
  s.set(239, -1);
  s.set(263, -1);
  s.set(503, -1);
  s.set(1261, -1);
  s.set(1501, -1);
  s.set(3783, -1);
  s.set(3927, -1);
  s.set(3933, -1);
  s.set(4565, -1);
  s.set(4685, 1);
  s.set(6233, 1);

  DenseRow<int,int> d(&s);
  BOOST_CHECK_EQUAL( d.first_nonzero_column(), 0 );
  for(int i = 0; i < 7260; ++i)
    if (6233 == i or 4685 == i or 0 == i)
      BOOST_CHECK_EQUAL(d.get(i), 1);
    else if (239  == i or
             263  == i or
             503  == i or
             1261 == i or
             1501 == i or
             3783 == i or
             3927 == i or
             3933 == i or
             4565 == i)
      BOOST_CHECK_EQUAL(d.get(i), -1);
    else 
      BOOST_CHECK_EQUAL(d.get(i), 0);
}


BOOST_AUTO_TEST_CASE( check_denserow_ge_with_sparse1 )
{
  DenseRow<int,int>* d = new DenseRow<int,int>(0,100);
  d->set(0, 1);
  d->set(90, 9);

  SparseRow<int,int>* s = new SparseRow<int,int>(0,100, 1);  

  // check elimination
  DenseRow<int,int>* dd = d->gaussian_elimination(s);
  BOOST_CHECK ( d != dd );
  BOOST_CHECK_EQUAL( dd->first_nonzero_column(), 90 );
  BOOST_CHECK_EQUAL( dd->get(dd->first_nonzero_column()), -9 );
  BOOST_CHECK_EQUAL( dd->size(), d->size() - 90 );
}


BOOST_AUTO_TEST_CASE( check_denserow_ge_with_sparse2 )
{
  // init d - only at odd places
  DenseRow<int,int>* d = new DenseRow<int,int>(0,100);
  d->set(0, 5);
  for(int i = 1; i <= 100; i += 2)
    d->set(i, 7);

  // init s - only at even places
  SparseRow<int,int> *s = new SparseRow<int,int>(0, 100, 2);
  for(int i = 2; i <= 100; i += 2)
    s->set(i, 3);

  // check elimination
  DenseRow<int,int>* dd = d->gaussian_elimination(s);
  BOOST_CHECK_EQUAL( dd->size(), d->size()-1 );
  BOOST_CHECK_EQUAL( dd->first_nonzero_column(), 1 );
  for (int i = 1; i <= 100; ++i)
    BOOST_CHECK_EQUAL(dd->get(i), (i%2 == 0)? 5*3 : -2*7);
}
