// Hey Emacs, this is -*- c++ -*-

#include <cstdlib>
#include <sstream>

rheinfall::Rank<val_t, coord_t>*
matrix_from_file(const std::string& filename)
{
  // automake stores the path to the source directory
  // into environment variable $srcdir
  std::ostringstream path;
  const char * const srcdir = getenv("srcdir");
  if (NULL != srcdir)
    path << getenv("srcdir") << "/";
  path << filename;

  std::ifstream input(path.str().c_str());
  if (input.fail()) {
    std::ostringstream msg;
    msg << "Cannot open file '" <<path.str()<< "' for reading - skipping." 
        << std::endl;
    throw std::runtime_error(msg.str());
  };
  
  // read matrix dimensions
  long rows, cols;
  char M;
  input >> std::skipws >> rows >> cols >> M;
  if ('M' != M) {
    std::ostringstream msg;
    msg << "File '" <<filename<< "' lacks SMS header - skipping." 
        << std::endl;
    throw std::runtime_error(msg.str());
  };

  rheinfall::Rank<val_t,coord_t>* m 
    = new rheinfall::Rank<val_t,coord_t>(cols, 2);

  m->read_triples(input, rows, cols, true, false);
  input.close();
  
  return m;
};


BOOST_AUTO_TEST_CASE( check_M05_D4 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,5-D4.sms");
  BOOST_CHECK_EQUAL( m->rank(), 210 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D5 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,5-D5.sms");
  BOOST_CHECK_EQUAL( m->rank(), 1902 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D6 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,5-D6.sms");
  BOOST_CHECK_EQUAL( m->rank(), 5358 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D7 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,5-D7.sms");
  BOOST_CHECK_EQUAL( m->rank(), 5916 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D8 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,5-D8.sms");
  BOOST_CHECK_EQUAL( m->rank(), 2239 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D5 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,6-D5.sms");
  BOOST_CHECK_EQUAL( m->rank(), 3024 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D6 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,6-D6.sms");
  BOOST_CHECK_EQUAL( m->rank(), 46776 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D7 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,6-D7.sms");
  BOOST_CHECK_EQUAL( m->rank(), 247704 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D8 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,6-D8.sms");
  BOOST_CHECK_EQUAL( m->rank(), 614586 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D11 )
{
  rheinfall::Rank<val_t,coord_t>* m = matrix_from_file("M0,6-D11.sms");
  BOOST_CHECK_EQUAL( m->rank(), 122879 );
  delete m;
}
