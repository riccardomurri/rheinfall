// Hey Emacs, this is -*- c++ -*-

rheinfall::Rheinfall<val_t, coord_t>*
matrix_from_file(const std::string& filename)
{
  std::ifstream input(filename.c_str());
  if (input.fail()) {
    std::ostringstream msg;
    msg << "Cannot open file '" <<filename<< "' for reading - skipping." 
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

  // possibly transpose matrix so that rows = min(rows, cols)
  // i.e., minimize the number of eliminations to perform
  bool transpose = false;
  if (rows > cols) {
    std::swap(rows, cols);
    transpose = true;
  };
  
  rheinfall::Rheinfall<val_t,coord_t>* m 
    = new rheinfall::Rheinfall<val_t,coord_t>(cols, 1);

  m->read(input, rows, cols, true, transpose);
  input.close();
  
  return m;
};


BOOST_AUTO_TEST_CASE( check_M05_D4 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M05-D4.txt");
  BOOST_CHECK_EQUAL( m->rank(), 210 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D5 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M05-D5.txt");
  BOOST_CHECK_EQUAL( m->rank(), 1902 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D6 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M05-D6.txt");
  BOOST_CHECK_EQUAL( m->rank(), 5358 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D7 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M05-D7.txt");
  BOOST_CHECK_EQUAL( m->rank(), 5916 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M05_D8 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M05-D8.txt");
  BOOST_CHECK_EQUAL( m->rank(), 2239 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D5 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M06-D5.txt");
  BOOST_CHECK_EQUAL( m->rank(), 3024 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D6 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M06-D6.txt");
  BOOST_CHECK_EQUAL( m->rank(), 46776 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D7 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M06-D7.txt");
  BOOST_CHECK_EQUAL( m->rank(), 247704 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D8 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M06-D8.txt");
  BOOST_CHECK_EQUAL( m->rank(), 614586 );
  delete m;
}


BOOST_AUTO_TEST_CASE( check_M06_D11 )
{
  rheinfall::Rheinfall<val_t,coord_t>* m = matrix_from_file("test/inputs/M06-D11.txt");
  BOOST_CHECK_EQUAL( m->rank(), 122879 );
  delete m;
}
