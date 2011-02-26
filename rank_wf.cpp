/**
 * @file   rank_wf.cpp
 *
 * Sample implementation of the ``waterfall'' algorithm for computing matrix rank.
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

// TODO:
// === first working version ===
//   - have read_*() return the number of nonzero entries
//   - check Processor phases to optimize communication/computation overlap


#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>

#include <sys/time.h>
#include <sys/resource.h>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
#include <boost/variant/get.hpp>
namespace mpi = boost::mpi;


enum { TAG_END=0, TAG_ROW_SPARSE=1, TAG_ROW_DENSE=2 };


typedef int coord_t;
typedef long long val_t;


/** Implements the `Waterfall` algorithm for distributed computation
    of matrix rank. */
class Waterfall {

public:
  /** Constructor, taking row and column dimension; uses the @c boost::mpi
      default communicator (@c MPI_COMM_WORLD ). */
  Waterfall(const coord_t nrows, const coord_t ncols);
  
  /** Read a matrix stream into the processors. Assumes that entries
      belonging to one row are not interleaved with entries from other
      rows; i.e., that it can consider reading row @i i complete as
      soon as it finds a row index @i j != i. Return number of nonzero
      entries _read_. */
  size_t read_noninterleaved(std::istream& input);

  /** Read a matrix stream into the processors. This should be run
      from one MPI rank only: it will then distribute the data to
      other MPI ranks. Return total number of nonzero entries read. */
  size_t read_and_distribute(std::istream& input) { throw "not implemented"; };

  /** Return rank of matrix after in-place destructive computation. */
  int rank();

  // forward declarations
protected:
  class Processor;
  

private:
  mpi::communicator comm;
  const int mpi_id; /**< MPI rank. */
  const int mpi_nprocs; /**< Total number of ranks in MPI communicator. */
  
  const coord_t nrows;
  const coord_t ncols;
  
  std::vector<Processor*> procs;   /**< Local processors. */

protected:
  // implement cyclic distribution of columns to MPI ranks
  bool is_local(const coord_t c) const; 
  Processor& local_owner(const coord_t c);
  int remote_owner(const coord_t c) const;

  /** A single processing element. */
  class Processor {

  public:
    Processor(Waterfall& parent_, 
              const coord_t starting_column_, const coord_t ending_column_,
              const double dense_threshold_ = 40.0) 
      : parent(parent_),
        starting_column(starting_column_), ending_column(ending_column_), rows(), 
        phase(running), inbox(), outbox(), 
        dense_threshold(dense_threshold_) 
    { };

  public:
    // forward declarations
    class SparseRow;
    class DenseRow;
    typedef SparseRow* sparse_row_ptr;
    typedef DenseRow* dense_row_ptr;
    typedef boost::variant< sparse_row_ptr,dense_row_ptr > row_ptr;
#ifndef NDEBUG
    static sparse_row_ptr __get_sparse_row_ptr(row_ptr r) {
      return boost::get<sparse_row_ptr,sparse_row_ptr,dense_row_ptr>(r); };
    static dense_row_ptr __get_dense_row_ptr(row_ptr r) { 
      return boost::get<dense_row_ptr,sparse_row_ptr,dense_row_ptr>(r); };
#endif

  protected:     
    Waterfall& parent;
    const coord_t starting_column;
    const coord_t ending_column;

    /** The block of rows owned by this processor */
    typedef std::vector<row_ptr> block;
    block rows;

    /** Processor state: it is `running` when Gaussian elimination
        operations are being carried out regularly; it turns @c ending
        when a @c TAG_END message is received; it is @c done when the
        final contribution to the rank has been computed and the @c
        step() method should be called no more. */
    enum { running, ending, done } phase;
    
    /** List of incoming messages from other processors; each message
        consists of a message tag (an int) and (optionally) a rows. */
    typedef std::list< row_ptr > inbox_t;
    inbox_t inbox;

    /** List of rows sent to other processors and the associated MPI request. 
        We need to keep track of these in order to free the resources when we're done. */
    typedef std::list< mpi::request > outbox_t;
    outbox_t outbox;

    /** A row will switch its storage to dense format when the percent of
        nonzero entries w.r.t. total length exceeds this value. */
    const double dense_threshold;

  public:
    class SparseRow 
    {
    public:
      SparseRow(const coord_t starting_column_, const coord_t ending_column_) 
        : starting_column(starting_column_), ending_column(ending_column_), storage() 
      { };
      /** Constructor initializing an invalid row; should call
          serialize immediately after!  */
      SparseRow() 
        : starting_column(-1), ending_column(-1), storage() 
      { };
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() const;
      /** Return number of allocated (i.e., non-zero) entries. */
      size_t size() const;
      /** Return a reference to the element stored at column @c col */
      val_t& operator[](const coord_t col);
      /** Return a copy of the element stored at column @c col */
      const val_t operator[](const coord_t col) const;
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. */
      sparse_row_ptr gaussian_elimination(sparse_row_ptr other);
      /** Perform Gaussian elimination: sum a multiple of this row and
          (a multiple of) row @c r so that the combination has its
          first nonzero entry at a column index strictly larger than
          the one of both rows. Return pointer to the combined row, which
          could possibly be this row if in-place update took place. 
          Elimination on a dense row returns a dense row. */
      dense_row_ptr gaussian_elimination(dense_row_ptr other);
    protected:
      coord_t starting_column; // would-be `const`: can only be modified by ctor and serialization
      coord_t ending_column; // would-be `const`: can only be modified by ctor and serialization
      val_t leading_term; // would-be `const`: can only be modified by ctor and serialization
      typedef std::vector< std::pair<coord_t,val_t> > storage_t;
      storage_t storage;
      friend class DenseRow;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive& ar, const unsigned int version);
    }; // class SparseRow

    class DenseRow
    {
    public:
      /** Constructor, copying entries from a @c SparseRow instance. */
      DenseRow(const sparse_row_ptr r);
      /** Construct a null row. */
      DenseRow(const coord_t starting_column_, const size_t ending_column_);
      /** Construct an invalid row; exists only for the purpose of
          calling @c serialize immediately after! */
      DenseRow(); // FIXME: make private?
      ~DenseRow();
      /** Return index of first column that has a nonzero entry. */
      coord_t first_nonzero_column() const;
      /** Return number of allocated entries. */
      size_t size() const;
      /** Return a reference to the element stored at column @c col */
      val_t& operator[](const coord_t col);
      /** Return a copy of the element stored at column @c col */
      const val_t operator[](const coord_t col) const;
      dense_row_ptr gaussian_elimination(const sparse_row_ptr other);
      dense_row_ptr gaussian_elimination(const dense_row_ptr other);
    protected:
      coord_t starting_column; 
      coord_t ending_column; 
      val_t leading_term; 
      /// Row length.  This is of course the number of elements in the @c storage array.
      size_t _size;
      /// Plain array of entries.
      val_t* storage; 
      friend class SparseRow;
      // serialization support (needed for MPI)
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive& ar, const unsigned int version);
      /** Adjust @c _size and @c starting_column to reflect the actual
          contents of the stored row. */
      size_t _resize();
    }; // class DenseRow

  public: 
    void recv_row(row_ptr& new_row);
    /** Send the row data pointed by @c row to the processor owning the
        block starting at @c col. Frees up @c row */
    void send_row(const coord_t column, row_ptr& row);
    class __sender : public boost::static_visitor<mpi::request> {
      mutable Waterfall& wf;
      const coord_t column;
    public:
      __sender(Waterfall& wf_, const coord_t column_);
      mpi::request operator() (sparse_row_ptr& r) const;
      mpi::request operator() (dense_row_ptr& r) const;
    };

    /** Make a pass over the block of rows and return either 0 or 1
        (depending whether elimination took place or not). */
    int step();
    class __delete_row 
      : public boost::static_visitor<void> {
    public:
      void operator() (sparse_row_ptr& r) const;
      void operator() (dense_row_ptr& r) const;
    }; // class __delet_row
    class __first_nonzero_column 
      : public boost::static_visitor<coord_t> {
    public:
      coord_t operator() (sparse_row_ptr& r) const;
      coord_t operator() (dense_row_ptr& r) const;
    }; // class __first_nonzero_column
    class __row_size 
      : public boost::static_visitor<size_t> {
    public:
      size_t operator() (const sparse_row_ptr& r) const;
      size_t operator() (const dense_row_ptr& r) const;
    }; // class __row_size



    void end_phase();

    bool is_done() const;

  protected:
    /** Send the termination signal to the processor owning the block
        starting at @c col. */
    void send_end(const coord_t column) const;
    /** Perform Gaussian elimination on row `second`, using the entry
        in the leading column of `first` as a pivot. */
    row_ptr gaussian_elimination(const row_ptr first, row_ptr second) const;
    class __ge_op
      : public boost::static_visitor<row_ptr>
    {
      const coord_t starting_column;
      const coord_t ending_column;
      const double dense_threshold;
    public:
      __ge_op(const coord_t starting_column_, const coord_t ending_column_, 
              const double dense_threshold_);
      row_ptr operator() (const sparse_row_ptr first, sparse_row_ptr second) const;
      row_ptr operator() (const dense_row_ptr first, sparse_row_ptr second) const;
      row_ptr operator() (const sparse_row_ptr first, dense_row_ptr second) const;
      row_ptr operator() (const dense_row_ptr first, dense_row_ptr second) const;
    }; // class __ge_op
  };
};


//
// implementation
//

/* Needed arithmetic operations */
static inline
val_t _gcd(long long const n, long long const m) 
{ 
  return m==0? n : _gcd(m, n%m);
};

static inline
val_t gcd(long long const n, long long const m) 
{ 
  return _gcd(std::abs(n), std::abs(m)); 
};


// local classes cannot be used as template argument; therefore we
// need to declare visitors for use with boost::variant at top-level
// :-(

inline void
Waterfall::Processor::__delete_row::operator() (sparse_row_ptr& r) const 
{
  delete r;
};

inline void
Waterfall::Processor::__delete_row::operator() (dense_row_ptr& r) const 
{
  delete r;
};


Waterfall::Processor::__ge_op::__ge_op(const coord_t starting_column_, 
                                       const coord_t ending_column_, 
                                       const double dense_threshold_)
  : starting_column(starting_column_), 
    ending_column(ending_column_),
    dense_threshold(dense_threshold_)
{ };

Waterfall::Processor::row_ptr 
Waterfall::Processor::__ge_op::operator() (const sparse_row_ptr first, 
                                           sparse_row_ptr second) const 
{
  // assume that `first` has already been checked for fill-in
  // in the main loop, so we won't attempt to convert it here;
  // just test `second` for too much fill-in
  const double fill_in = 100.0 * second->size() / (ending_column - starting_column + 1);
  if (fill_in > dense_threshold) {
    // FIXME: could merge ctor+elimination in one funcall
    dense_row_ptr dense_second(new DenseRow(second));
    delete second;
    return first->gaussian_elimination(dense_second);
  }
  else {
    return first->gaussian_elimination(second);
  }; // if (fill_in > dense_threshold)
};

Waterfall::Processor::row_ptr 
Waterfall::Processor::__ge_op::operator() (const dense_row_ptr first, 
                                           sparse_row_ptr second) const 
{
  return first->gaussian_elimination(second);
};
  
Waterfall::Processor::row_ptr 
Waterfall::Processor::__ge_op::operator() (const sparse_row_ptr first, 
                                           dense_row_ptr second) const 
{
  return first->gaussian_elimination(second);
};

Waterfall::Processor::row_ptr 
Waterfall::Processor::__ge_op::operator() (const dense_row_ptr first, 
                                           dense_row_ptr second) const 
{
  return first->gaussian_elimination(second);
};


inline coord_t 
Waterfall::Processor::__first_nonzero_column::operator() (sparse_row_ptr& r) const 
{
  return r->first_nonzero_column();
};

inline coord_t 
Waterfall::Processor::__first_nonzero_column::operator() (dense_row_ptr& r) const 
{
  return r->first_nonzero_column();
};


inline size_t 
Waterfall::Processor::__row_size::operator() (const sparse_row_ptr& r) const 
{
  return r->size();
};

inline size_t 
Waterfall::Processor::__row_size::operator() (const dense_row_ptr& r) const 
{
  return r->size();
};


Waterfall::Processor::__sender::__sender(Waterfall& wf_, const coord_t column_) 
  : wf(wf_), column(column_) 
{ };

inline mpi::request 
Waterfall::Processor::__sender::operator() (sparse_row_ptr& r) const
{
  return wf.comm.isend(wf.remote_owner(column), TAG_ROW_SPARSE, r);
};

inline mpi::request 
Waterfall::Processor::__sender::operator() (dense_row_ptr& r) const
{
  return wf.comm.isend(wf.remote_owner(column), TAG_ROW_DENSE, r);
};



// ------- Waterfall -------

Waterfall::Waterfall(const coord_t nrows, const coord_t ncols)
  : comm(mpi::communicator()), mpi_id(comm.rank()), mpi_nprocs(comm.size()),
    nrows(nrows), ncols(ncols),
    procs()
{  
  // setup the array of processors
  procs.reserve(1 + ncols / mpi_nprocs);
  for (int n = 0; n < 1 + (ncols / mpi_nprocs); ++n)
    procs.push_back(new Processor(*this, mpi_id + n * mpi_nprocs, ncols-1));
};


size_t
Waterfall::read_noninterleaved(std::istream& input) 
{
    // count non-zero items read
    size_t nnz = 0;

    // only read rows whose leading column falls in our domain
    std::vector< std::pair<coord_t,val_t> > row;

    coord_t current_row_index = -1;
    coord_t i, j;
    val_t value;
    while (not input.eof()) {
      input >> i >> j >> value;
      if (0 == i and 0 == j and 0 == value)
        break; // end of matrix stream
      ++nnz;
      // SMS indices are 1-based
      --i;
      --j;
      // new row?
      if (i != current_row_index) { // yes, commit old row
        // find starting column
        coord_t starting_column = ncols;
        for (std::vector< std::pair<coord_t,val_t> >::const_iterator it = row.begin();
             it != row.end();
             ++it)
          if (it->first < starting_column)
            starting_column = it->first;
        // commit row
        if (starting_column < ncols // otherwise, row is empty -- ignore it
            and is_local(starting_column)) {
          Processor::SparseRow* r = new Processor::SparseRow(starting_column, ncols-1);
          for (std::vector< std::pair<coord_t,val_t> >::const_iterator it = row.begin();
               it != row.end();
               ++it)
            {
              (*r)[it->first] = it->second;
            }; 
          boost::variant<Processor::SparseRow*,Processor::DenseRow*> r_(r);
          local_owner(starting_column).recv_row(r_);
        }; // if is_local(starting_column)
        row.clear();
        current_row_index = i;
      }; // if i != current_row_index
      row.push_back(std::make_pair(j, value));
    }; // while (! eof)
    
    return nnz;
};


int 
Waterfall::rank() 
{
  // kickstart termination signal
  if (is_local(0))
    local_owner(0).end_phase();
  // collect (partial) ranks
  std::vector<int> r(procs.size(), 0);
  size_t n0 = 0;
  while(true) {
    // n0 incremented each time a processor goes into `done` state 
    while (procs[n0]->is_done())
      ++n0;
    if (procs.size()-1 == n0)
      break;
    // step all processors
    for (size_t n = n0; n < procs.size(); ++n) {
      r[n] = procs[n]->step();
    };
    // manage arrival of messages
    while (boost::optional<mpi::status> status = comm.iprobe()) { 
      switch(status->tag()) {
      case TAG_ROW_SPARSE: {
        Processor::sparse_row_ptr new_row = new Processor::SparseRow();
        comm.recv(status->source(), status->tag(), new_row);
        assert(is_local(new_row->first_nonzero_column()));
        Processor::row_ptr new_row_ptr(new_row);
        local_owner(new_row->first_nonzero_column()).recv_row(new_row_ptr);
        break;
      };
      case TAG_ROW_DENSE: {
        Processor::dense_row_ptr new_row = new Processor::DenseRow();
        comm.recv(status->source(), status->tag(), 
                  *boost::get<Processor::DenseRow*>(new_row));
        assert(is_local(new_row->first_nonzero_column()));
        Processor::row_ptr new_row_ptr(new_row);
        local_owner(new_row->first_nonzero_column()).recv_row(new_row_ptr);
        break;
      };
      case TAG_END: {
        // "end" message received; no new rows will be coming.
        // But some other rows could have arrived or could
        // already be in the `inbox`, so we need to make another
        // pass anyway.  All this boils down to: set a flag,
        // make another iteration, and end the loop next time.
        coord_t column = -1;
        comm.recv(status->source(), status->tag(), column);
        assert(column >= 0 and is_local(column));
        local_owner(column).end_phase();
      };
      }; // switch(status->tag())
    };
  };

  // the partial rank is computed as the sum of all ranks computed
  // by local processors
  int partial_rank = 0;
  for (size_t n = 0; n < r.size(); ++n)
    partial_rank += r[n];
  
  // wait until all processors have done running
  comm.barrier();
  
  // collect the partial ranks for all processes
  int rank = 0;
  mpi::all_reduce(comm, partial_rank, rank, std::plus<int>());
  return rank;
};


inline bool 
Waterfall::is_local(const coord_t c) const 
{ 
  return (mpi_id == c % mpi_nprocs); 
};


inline Waterfall::Processor& 
Waterfall::local_owner(const coord_t c)
{ 
  return *(procs[c / mpi_nprocs]); 
};


inline int 
Waterfall::remote_owner(const coord_t c) const
{ 
  return (c % mpi_nprocs); 
};


// ------- Processor -------

inline void 
Waterfall::Processor::recv_row(row_ptr& new_row) 
{
  // FIXME: requires locking for MT operation
  inbox.push_back(new_row);
  // DEBUG
  sparse_row_ptr r = __get_sparse_row_ptr(new_row);
  std::cerr << "DEBUG: Processor " <<this 
            << " received row"
            << " of size " << r->size()
            << " starting at column " << r->first_nonzero_column()
            << " (sparse_row_ptr " <<r<< ")"
            << std::endl;
};


inline void 
Waterfall::Processor::send_row(const coord_t column, row_ptr& row) 
{
  if (parent.is_local(column))
    parent.local_owner(column).recv_row(row);
  else { // ship to remote process
    const __sender sender(parent, column);
    mpi::request req = boost::apply_visitor(sender, row);
    outbox.push_back(req);
  };
};


int 
Waterfall::Processor::step() 
{
  assert(not is_done()); // cannot be called once finished
  std::cerr << "DEBUG: stepping Processor " << this << std::endl;

  // delete row with boost::variant
  class _delete_row {
  public:
    void operator() (row_ptr& r) {
      return boost::apply_visitor(__delete_row(), r);
    };
  } delete_row;

  // get first nonzero column with boost::variant
  class _first_nonzero_column {
  public:
    coord_t operator() (row_ptr& r) {
      return boost::apply_visitor(__first_nonzero_column(), r);
    };
  } first_nonzero_column;

  // get row size with boost::variant
  class _row_size {
  public:
    size_t operator() (row_ptr& r) const {
      return boost::apply_visitor(__row_size(), r);
    };
  } row_size; 

  // receive new rows
  for (inbox_t::iterator it = inbox.begin(); it != inbox.end(); ++it) {
    row_ptr new_row = *it;
    rows.push_back(new_row);
    if (row_size(new_row) < row_size(rows.front()))
      // less nonzero entries: we have a new leading row
      std::swap(rows.front(), rows.back());
  };
  inbox.clear();
          
  if (rows.size() > 0) {
    row_ptr first = rows[0];
    for (block::iterator second = rows.begin()+1; second < rows.end(); ++second) {
      // perform elimination
      row_ptr new_row = gaussian_elimination(first, *second);
          
      // ship reduced rows to other processors
      if (row_size(new_row) > 0) {
        send_row(first_nonzero_column(new_row), new_row);
      }
      delete_row(new_row);
    };
    // just keep the "first" row
    rows.clear();
    rows.push_back(first);
  };

  if (running == phase) {
    if (outbox.size() > 0) {
      // check if some test messages have arrived
      // and free corresponding resources
      outbox.erase(mpi::test_some(outbox.begin(), outbox.end()), outbox.end());
    };
  }
  else { // `phase == ending`: end message already received
    // pass end message along
    if (starting_column < ending_column)
      send_end(starting_column + 1);
        
    if (outbox.size() > 0) {
      // wait untill all sent messages have arrived
      mpi::wait_all(outbox.begin(), outbox.end());
    };
        
    // all done
    phase = done;
  };
      
  // exit with 0 only if we never processed any row
  assert(rows.size() <= 1);
  return rows.size();
};


inline void 
Waterfall::Processor::end_phase() 
{ 
  phase = ending; 
};


inline bool 
Waterfall::Processor::is_done() const 
{ 
  return (done == phase); 
};


inline void 
Waterfall::Processor::send_end(const coord_t column) const
{
  if (parent.is_local(column))
    parent.local_owner(column).end_phase();
  else { // ship to remote process
    parent.comm.send(parent.remote_owner(column), TAG_END, column);
  };
}; // send_end


inline Waterfall::Processor::row_ptr 
Waterfall::Processor::gaussian_elimination(const row_ptr first, row_ptr second) const
{
  const __ge_op ge_op(starting_column, ending_column, dense_threshold);
  return boost::apply_visitor(ge_op, first, second);
}; // gaussian_elimination


// ------- SparseRow -------

inline coord_t 
Waterfall::Processor::SparseRow::first_nonzero_column() const 
{
  assert(0 != leading_term);
  return starting_column;
};

inline size_t 
Waterfall::Processor::SparseRow::size() const 
{ 
  return storage.size(); 
};


Waterfall::Processor::sparse_row_ptr 
Waterfall::Processor::SparseRow::gaussian_elimination(sparse_row_ptr other) 
{
  assert(this->starting_column == other->starting_column);

  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term, other->leading_term);
  val_t a = - (other->leading_term / GCD);
  val_t b = leading_term / GCD;

  sparse_row_ptr result = NULL; // XXX: use boost::optional<...> instead?

  storage_t::const_iterator this_i = this->storage.begin(); 
  storage_t::const_iterator other_i = other->storage.begin();
  // loop while one of the two indexes is still valid
  while(this_i != this->storage.end() or other_i != other->storage.end()) {
    coord_t this_col;
    if (this_i != this->storage.end()) {
      this_col = this_i->first;
      assert (this_col > this->starting_column);
    }
    else // this_i reached end of vector, use out-of-range value
      this_col = this->ending_column + 1;

    coord_t other_col;
    if (other_i != other->storage.end()) {
      other_col = other_i->first;
      assert (other_col > other->starting_column);
    }
    else 
      other_col = other->ending_column + 1;

    bool nonzero_found = false;
    coord_t coord;
    val_t entry;
    if (other_col < this_col) {
      entry = b * other_i->second;
      if (0 != entry) {
        coord = other_col;
        nonzero_found = true;
      };
      ++other_i;
    }
    else if (other_col == this_col) {
      entry = a * this_i->second + b * other_i->second;
      if (0 != entry) {
        coord = this_col;
        nonzero_found = true;
      };
      ++this_i;
      ++other_i;
    }
    else if (other_col > this_col) {
      entry = a * this_i->second;
      if (0 != entry) {
        coord = this_col;
        nonzero_found = true;
      };
      ++this_i;
    };
    
    if (nonzero_found) {
      // XXX: the following code makes assumptions about how the
      // storage is organised within `SparseRow`; it duplicates code
      // from SparseRow::operator[] for efficiency
      if (NULL == result) {
        // allocate new SparseRow
        result = new SparseRow(coord, ending_column);
        result->leading_term = entry;
        result->storage.reserve(storage.size() + other->storage.size());
      }
      else {
        // FIXME: change storage order so that this is push_back()!
        result->storage.insert(result->storage.begin(), std::make_pair(coord, entry));
      }; 
    }; 
  }; // while(this_i < storage.rend() ...)
  delete other; // release old storage
  if (NULL == result) // row full of zeroes
    return new SparseRow(ending_column, ending_column);
  else
    return result;
}; // sparse_row_ptr gaussian_elimination(sparse_row_ptr r)


Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::SparseRow::gaussian_elimination(dense_row_ptr other) 
{
  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term, other->leading_term);
  val_t a = - (other->leading_term / GCD);
  val_t b = leading_term / GCD;

  for (size_t j = 0; j < other->size(); ++j) 
    other->storage[j] *= b;
  for (storage_t::const_iterator it = storage.begin();
       it != storage.end();
       ++it)
    {
      assert(it->first > starting_column and it->first <= ending_column);
      other->storage[it->first] += a*it->second;
    };

  other->_resize(); // content update done, adjust size and starting column
  return other;
}; // dense_row_ptr gaussian_elimination(dense_row_ptr other)


inline val_t& 
Waterfall::Processor::SparseRow::operator[](const coord_t col) 
{
  assert(col >= starting_column and col <= ending_column);
  if (col == starting_column) 
    return leading_term;
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < storage.size() and storage[jj].first < col)
    ++jj;
  if (storage.size() == jj) {
    // end of list reached, `col` is larger than any index in this row
    storage.push_back(std::make_pair<size_t,val_t>(col,0));
    return storage.back().second;
  }
  else if (col == storage[jj].first) 
    return storage[jj].second;
  else { // storage[jj].first > j > storage[jj+1].first
    // insert new pair before `jj+1`
    storage.insert(storage.begin()+jj+1, 
                   std::make_pair<size_t,val_t>(col,0));
    return storage[jj+1].second;
  };
}; // SparseRow::operator[](...)


inline const val_t
Waterfall::Processor::SparseRow::operator[](const coord_t col) const
{
  assert(col >= starting_column and col <= ending_column);
  if (col == starting_column) 
    return leading_term;
  if (storage.size() == 0)
    return 0;
  // else, fast-forward to place where element is/would be stored
  size_t jj = 0;
  while (jj < storage.size() and storage[jj].first < col)
    ++jj;
  if (storage.size() == jj) {
    // end of list reached, `col` is larger than any index in this row
    return 0;
  }
  else if (col == storage[jj].first) 
    return storage[jj].second;
  else { // storage[jj].first > j > storage[jj+1].first
    return 0;
  };
}; // SparseRow::operator[](...) const 


template<class Archive>
void 
Waterfall::Processor::SparseRow::serialize(Archive& ar, const unsigned int version) 
{
  assert(version == 0);
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  ar & starting_column & ending_column & leading_term & storage;
  assert(starting_column >= 0);
  assert(ending_column >= 0);
}; // SparseRow::serialize(...)


// ------- DenseRow -------

Waterfall::Processor::DenseRow::DenseRow(const sparse_row_ptr r) 
  : starting_column(r->starting_column), 
    ending_column(r->ending_column),
    _size(ending_column - starting_column + 1),
    storage(new val_t[_size]) 
{ 
  assert(std::distance(r->storage.begin(), r->storage.end()) <= _size);
  std::fill_n(storage, _size, 0);
  for (SparseRow::storage_t::const_iterator it = r->storage.begin();
       it != r->storage.end();
       ++it) 
    {
      storage[it->first] = it->second;
    };
};


inline
Waterfall::Processor::DenseRow::DenseRow(const coord_t starting_column_, 
                                         const size_t ending_column_) 
  : starting_column(starting_column_), 
    ending_column(ending_column_),
    _size(ending_column - starting_column + 1),
    storage(new val_t[_size]) 
{ 
  std::fill_n(storage, _size, 0);
};


inline
Waterfall::Processor::DenseRow::DenseRow() 
  : starting_column(-1), 
    ending_column(-1),
    _size(0),
    storage(NULL) 
{ };


inline
Waterfall::Processor::DenseRow::~DenseRow() 
{ 
  delete [] storage; 
};


inline coord_t 
Waterfall::Processor::DenseRow::first_nonzero_column() const 
{
  assert (0 != leading_term);
  return starting_column;
};


inline Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::DenseRow::gaussian_elimination(const sparse_row_ptr other) 
{
  // convert `other` to dense storage upfront: adding the
  // non-zero entries from `this` would made it pretty dense
  // anyway
  dense_row_ptr dense_other(new DenseRow(other));
  delete other;
  return this->gaussian_elimination(dense_other);
}; // dense_row_ptr gaussian_elimination(sparse_row_ptr other)


inline Waterfall::Processor::dense_row_ptr 
Waterfall::Processor::DenseRow::gaussian_elimination(const dense_row_ptr other) 
{
  assert(this->size() == other->size());
  // compute:
  //   `a`: multiplier for `this` row
  //   `b`: multiplier for `other` row
  const val_t GCD = gcd(leading_term, other->leading_term);
  val_t a = - (other->leading_term / GCD);
  val_t b = leading_term / GCD;

  for (size_t j = 0; j < size(); ++j) {
    // XXX: is it faster to allocate new storage and fill it with `a*x+b*y`?
    other->storage[j] *= a;
    other->storage[j] += b * storage[j];
  };
  other->_resize(); // update done, adjust size and string column
  return other;
}; // dense_row_ptr gaussian_elimination(dense_row_ptr other)


inline val_t& 
Waterfall::Processor::DenseRow::operator[](const coord_t col) 
{
#ifndef NDEBUG
  if (col < starting_column or col > ending_column)
    throw std::out_of_range("DenseRow::operator[] was passed a column index out of valid range");
#endif
  if (col == starting_column)
    return leading_term;
  else
    return storage[_size - col - 1];
};


inline const val_t
Waterfall::Processor::DenseRow::operator[](const coord_t col) const
{
#ifndef NDEBUG
  if (col < starting_column or col > ending_column)
    throw std::out_of_range("DenseRow::operator[] was passed a column index out of valid range");
#endif
  if (col == starting_column)
    return leading_term;
  else
    return storage[_size - col - 1];
};


inline size_t
Waterfall::Processor::DenseRow::_resize()
{
  // compute new starting column
  coord_t new_starting_column = ending_column;
  assert(_size > 0);
  for (size_t j = 0; j < _size; ++j)
    if ((*this)[j] != 0)
      new_starting_column = starting_column + j;

  _size -= new_starting_column - starting_column;
  starting_column = new_starting_column;
  return _size;
};


template<class Archive>
void 
Waterfall::Processor::DenseRow::serialize(Archive& ar, const unsigned int version) 
{
  assert(version == 0);
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  ar & starting_column & ending_column & leading_term & _size;
  assert(starting_column >= 0);
  assert(ending_column >= 0);
  assert(_size >= 0);
  if (NULL == storage)
    storage = new val_t[_size];
  for (size_t j = 0; j < _size; ++j) {
    ar & storage[j];
  };
};


inline size_t 
Waterfall::Processor::DenseRow::size() const 
{ 
  return _size; 
};


//
// main
//

int
main(int argc, char** argv)
{
  if (1 == argc 
      or std::string("--help") == argv[1] 
      or std::string("-h") == argv[1]) 
    {
      std::cout << "Usage: " << argv[0] << " matrix_file.sms" << std::endl;
      return 1;
    }

  mpi::environment env(argc, argv);
  mpi::communicator world;
  int myid = world.rank();

  for (int i = 1; i < argc; ++i)
    {
      std::ifstream input(argv[i]);
      if (input.fail()) {
        std::cerr << "Cannot open file '" <<argv[i]<< "' for reading - skipping." 
                  << std::endl;
        continue;
      };

      if (0 == myid)
        std::cout << argv[0] << " file:"<<argv[i];

      // read matrix dimensions
      size_t rows, cols;
      char M;
      input >> std::skipws >> rows >> cols >> M;
      assert ('M' == M);
      if (0 == myid) {
        std::cout << " rows:" << rows;
        std::cout << " cols:" << cols;
      };

      Waterfall w(rows, cols);
      size_t nnz = w.read_noninterleaved(input);
      input.close();
      if (0 == myid)
        std::cout << " nonzero:" << nnz;

      struct rusage ru;
      getrusage(RUSAGE_SELF, &ru);
      struct timeval t0; memcpy(&t0, &ru.ru_utime, sizeof(struct timeval));

      if (0 == myid)
        std::cout << " rank:" << w.rank();

      getrusage(RUSAGE_SELF, &ru);
      struct timeval t1; memcpy(&t1, &ru.ru_utime, sizeof(struct timeval));
      struct timeval tdelta; timersub(&t1, &t0, &tdelta);
      double elapsed = tdelta.tv_sec + (tdelta.tv_usec / 1000000.0);
      if (0 == myid)
        std::cout << " time:" << std::fixed << std::setprecision(2) << elapsed << std::endl;
    }

  return 0;
}


